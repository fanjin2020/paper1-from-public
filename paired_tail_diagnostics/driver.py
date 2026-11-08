"""Targeted, fresh-seed tail and weight-source controls; old experiments unchanged."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime
import gzip
import math
import os
from pathlib import Path
import shutil
import signal
import statistics
import time
import uuid
import zipfile
from paired_count_calibration import driver as b
from . import VERSION

ROOT=b.ROOT
ARMS=('estimated_shared','estimated_design','independent_shared','independent_design','fixed_unit','gaussian_unit')
KEYS={(a,t) for a in ARMS for t in ('model_based','HC3')}
FILES=tuple(dict.fromkeys(b.FILES+('paper1_paired_tail.py','run_paired_tail_diagnostics.sh',
    'paired_tail_diagnostics/__init__.py','paired_tail_diagnostics/driver.py',
    'paired_tail_diagnostics/engine.R','paired_tail_diagnostics/worker.R','paired_tail_diagnostics/README.md')))
LABELS={'p1e6','p1e5','bh_first','p5e5','p1e4','p1e3','p01','p05'}


def scenarios():
    return [s for s in b.scenarios('full') if s['id'] in {'s01','s09','s11','s12','s15','s23'}]


def validate(r,p,s,rep):
    expected=p['settings']['seed']+s['index']*100000+rep
    if not isinstance(r,dict) or (r.get('run_id'),r.get('scenario'),r.get('replicate'),r.get('seed'))!=(p['run_id'],s['id'],rep,expected):
        raise ValueError('Result identity mismatch')
    if r.get('independent_seed')!=expected+50000000 or r.get('gaussian_seed')!=expected+100000000:
        raise ValueError('Independent random streams differ')
    rows=r.get('metrics',[])
    if len(rows)!=12 or {(a['weight_mode'],a['test']) for a in rows}!=KEYS:
        raise ValueError('Missing or duplicate arm')
    for a in rows:
        g,k,n=a['genes'],a['retained'],a['available']
        if any(type(x) is not int for x in (g,k,n)) or g!=p['settings']['genes'] or not 0<n<=k<=g:
            raise ValueError('Invalid availability')
        if a['filtered']!=g-k or a['unavailable_retained']!=k-n or a['paired_truth']!=0:
            raise ValueError('Invalid counts or truth')
        if a['response']!=('independent_gaussian' if a['weight_mode']=='gaussian_unit' else 'observed_logCPM'):
            raise ValueError('Wrong response identity')
        for key in ('coverage95','fpr_050','min_p'):
            if not b.finite(a[key]) or not 0<=a[key]<=1:raise ValueError('Invalid '+key)
        if not b.finite(a['lmfit_max_difference']) or a['lmfit_max_difference']>=1e-7:
            raise ValueError('lmFit/paired contrast mismatch')
        for family in ('planned','retained'):
            for level in ('010','050','100'):
                count=a[f'discoveries_{family}_{level}']
                if type(count) is not int or not 0<=count<=n or a[f'fdp_{family}_{level}']!=int(count>0):
                    raise ValueError('Invalid global-null FDP')
        if len(a['tails'])!=8 or {t['label'] for t in a['tails']}!=LABELS:raise ValueError('Missing tail thresholds')
        cutoffs=dict(p1e6=1e-6,p1e5=1e-5,bh_first=.05/g,p5e5=5e-5,p1e4=1e-4,p1e3=1e-3,p01=.01,p05=.05)
        previous=-1
        for t in sorted(a['tails'],key=lambda t:t['cutoff']):
            if abs(t['cutoff']-cutoffs[t['label']])>1e-15 or type(t['count']) is not int or not previous<=t['count']<=n:
                raise ValueError('Invalid tail counts')
            previous=t['count']
        if abs(next(t['count'] for t in a['tails'] if t['label']=='p05')-a['fpr_050']*n)>1e-8:
            raise ValueError('Strict tail/FPR mismatch')
        if sum(a['p_histogram'])!=n or abs(a['p_histogram'][0]-a['fpr_050']*n)>1e-8:
            raise ValueError('Histogram boundary mismatch')
        if sum(v['available'] for v in a['strata'])!=n or sum(v['retained'] for v in a['strata'])!=k:
            raise ValueError('Incomplete expression strata')
    return r


def valid(p,s,rep):
    d=b.folder(p,s,rep);receipt=b.read(d/'receipt.json')
    try:
        if not isinstance(receipt,dict) or receipt.get('run_id')!=p['run_id'] or set(receipt['files'])!={'result.json','events.csv.gz'}:return False
        for name,info in receipt['files'].items():
            if (d/name).stat().st_size!=info['bytes'] or b.sha(d/name)!=info['sha256']:return False
        validate(b.read(d/'result.json'),p,s,rep)
        return True
    except (OSError,ValueError,KeyError,TypeError):return False


def receipt(p,s,rep):
    d=b.folder(p,s,rep);r=validate(b.read(d/'result.json'),p,s,rep)
    with gzip.open(d/'events.csv.gz','rt',encoding='utf-8',newline='') as f:
        events=list(csv.DictReader(f))
    counts={key:0 for key in KEYS}
    seen=set()
    for a in events:
        key=(a['arm'],a['test']);identity=(*key,a['gene_id'])
        if key not in KEYS or identity in seen or a['scenario']!=s['id'] or int(a['replicate'])!=rep:
            raise ValueError('Bad event rows')
        if not 0<=float(a['p'])<=1:raise ValueError('Bad event p')
        if not 0<=float(a['q_retained'])<=float(a['q_planned'])+1e-12<=1+1e-12:
            raise ValueError('Bad event q values')
        seen.add(identity);counts[key]+=1
    for a in r['metrics']:
        if counts[(a['weight_mode'],a['test'])]!=a['event_candidates']-a['events_truncated']:
            raise ValueError('Truncated or missing events file')
    b.atomic_json(d/'receipt.json',dict(run_id=p['run_id'],files={n:dict(bytes=(d/n).stat().st_size,sha256=b.sha(d/n)) for n in ('result.json','events.csv.gz')}))


def progress(p):
    pending=[];done=0
    for i,s in enumerate(p['scenarios'],1):
        missing=[r for r in range(1,p['settings']['reps']+1) if not valid(p,s,r)]
        done+=p['settings']['reps']-len(missing);pending.append((i,s,missing))
    return done,pending


def report(p):
    out=Path(p['output']);done,_=progress(p)
    if done!=len(p['scenarios'])*p['settings']['reps']:raise ValueError('Cannot summarize incomplete run')
    groups={};flat=[];warnings=[];tails=[];strata=[];pairs=[];summary=[]
    event_tmp=out/'tail_events.csv.gz.tmp'
    with gzip.open(event_tmp,'wt',encoding='utf-8',newline='') as ef:
        writer=None
        for s in p['scenarios']:
            for rep in range(1,p['settings']['reps']+1):
                r=b.read(b.folder(p,s,rep)/'result.json')
                warnings.extend(dict(scenario=s['id'],replicate=rep,message=w) for w in r['warnings'])
                for a in r['metrics']:
                    key=(s['id'],a['weight_mode'],a['test']);groups.setdefault(key,[]).append(a)
                    row={k:v for k,v in a.items() if k not in {'tails','strata','p_histogram'}}
                    row.update(scenario=s['id'],replicate=rep,seed=r['seed'])
                    flat.append(row)
                with gzip.open(b.folder(p,s,rep)/'events.csv.gz','rt',encoding='utf-8',newline='') as f:
                    reader=csv.DictReader(f)
                    if writer is None:writer=csv.DictWriter(ef,reader.fieldnames);writer.writeheader()
                    writer.writerows(reader)
    os.replace(event_tmp,out/'tail_events.csv.gz')
    for s in p['scenarios']:
        for arm,test in sorted(KEYS):
            a=groups[(s['id'],arm,test)]
            row=dict(scenario=s['id'],n=s['n'],rho=s['rho'],depth=s['depth'],exposure_fraction=s['exposure_fraction'],arm=arm,test=test,datasets=len(a))
            for k in ('coverage95','fpr_050','mean_delta','mean_se','retained','available','events_truncated'):
                row[k],row[k+'_mcse']=b.mean_mc([v[k] for v in a])
            for family in ('planned','retained'):
                k=f'fdp_{family}_050';events=sum(v[k] for v in a)
                row[k]=events/len(a);row[k+'_ci_low'],row[k+'_ci_high']=b.wilson(events,len(a))
            bonf=sum(next(t['count'] for t in v['tails'] if t['label']=='bh_first')>0 for v in a)
            row['bonferroni_fwer']=bonf/len(a)
            row['bonferroni_ci_low'],row['bonferroni_ci_high']=b.wilson(bonf,len(a))
            summary.append(row)
            for label in sorted(LABELS):
                ts=[next(t for t in v['tails'] if t['label']==label) for v in a]
                rates=[t['count']/v['available'] for t,v in zip(ts,a)]
                mean,se=b.mean_mc(rates);cut=ts[0]['cutoff']
                any_count=sum(t['count']>0 for t in ts)
                any_low,any_high=b.wilson(any_count,len(a))
                tails.append(dict(scenario=s['id'],arm=arm,test=test,label=label,cutoff=cut,
                    exceedances=sum(t['count'] for t in ts),available=sum(v['available'] for v in a),
                    mean_rate=mean,rate_mcse=se,ratio_to_uniform=mean/cut,ratio_mcse=se/cut if se is not None else None,
                    datasets_with_any=any_count,any_probability=any_count/len(a),any_ci_low=any_low,any_ci_high=any_high))
            for label in ('mean_lt10','mean10_lt50','mean_ge50'):
                st=[next(t for t in v['strata'] if t['stratum']==label) for v in a]
                row=dict(scenario=s['id'],arm=arm,test=test,stratum=label)
                for k in ('retained','available','raw05','tail1e4','discoveries_planned','discoveries_retained'):
                    row[k]=sum(v[k] for v in st)
                for k in ('raw05','tail1e4'):
                    rates=[v[k]/v['available'] for v in st if v['available']>0]
                    row[k+'_mean_rate'],row[k+'_mcse']=b.mean_mc(rates) if rates else (None,None)
                strata.append(row)
        # Only compare arms sharing the observed response. Gaussian is a separate control.
        comparisons=[(a,t,'estimated_shared',t) for a in ARMS[1:-1] for t in ('model_based','HC3')]
        comparisons += [(a,'HC3',a,'model_based') for a in ARMS]
        for arm,test,ref,rt in comparisons:
            a=groups[(s['id'],arm,test)];z=groups[(s['id'],ref,rt)]
            for metric in ('fdp_planned_050','fdp_retained_050','fpr_050','coverage95'):
                delta,se=b.mean_mc([x[metric]-y[metric] for x,y in zip(a,z)])
                pairs.append(dict(scenario=s['id'],arm=arm,test=test,reference_arm=ref,reference_test=rt,
                    metric=metric,mean_difference=delta,mcse=se,
                    low=delta-1.96*se if se is not None else None,high=delta+1.96*se if se is not None else None))
    for name,rows in [('summary.csv',summary),('tail_summary.csv',tails),('expression_strata.csv',strata),('paired_comparisons.csv',pairs)]:
        b.atomic_csv(out/name,rows)
    with gzip.open(out/'replicate_metrics.csv.gz.tmp','wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(flat[0]));w.writeheader();w.writerows(flat)
    os.replace(out/'replicate_metrics.csv.gz.tmp',out/'replicate_metrics.csv.gz')
    b.atomic_json(out/'warnings.json',warnings)
    b.atomic_json(out/'completion.json',dict(run_id=p['run_id'],status='complete',datasets=done,
        method_rows=len(flat),warnings=len(warnings),scope='diagnostic controls, not a calibrated replacement test'))
    names=['protocol.json','environment.json','completion.json','summary.csv','tail_summary.csv','expression_strata.csv',
        'paired_comparisons.csv','replicate_metrics.csv.gz','tail_events.csv.gz','warnings.json']
    b.atomic_json(out/'review_manifest.json',{n:b.sha(out/n) for n in names})
    with zipfile.ZipFile(out/'review_tail.zip.tmp','w',compression=zipfile.ZIP_DEFLATED) as z:
        for n in names+['review_manifest.json']:z.write(out/n,n)
        for n in FILES:z.write(ROOT/n,'source/'+n)
    os.replace(out/'review_tail.zip.tmp',out/'review_tail.zip')
    print('COMPLETE',done,'datasets\nREVIEW',out/'review_tail.zip',flush=True)


def main(argv=None):
    def terminate(signum,frame):raise KeyboardInterrupt
    signal.signal(signal.SIGTERM,terminate)
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path);ap.add_argument('--workers',type=int,default=36)
    ap.add_argument('--profile',choices=('smoke','pilot','full'));ap.add_argument('--reps',type=int);ap.add_argument('--genes',type=int)
    ap.add_argument('--seed',type=int);ap.add_argument('--resume',action='store_true');ap.add_argument('--status',action='store_true')
    ap.add_argument('--collect',action='store_true');ap.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
    ap.add_argument('--r-library',action='append',default=[]);ap.add_argument('--allow-version-drift',action='store_true')
    args=ap.parse_args(argv)
    if args.workers<1:ap.error('Workers must be positive')
    if (args.resume or args.status or args.collect) and args.output is None:ap.error('Need actual --output directory')
    out=(args.output or ROOT.parent/('p1_tail_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).expanduser().resolve()
    if out==ROOT or out.is_relative_to(ROOT) or ROOT.is_relative_to(out):ap.error('Output overlaps code')
    p=b.read(out/'protocol.json')
    if p and p.get('experiment')!='paired_tail_diagnostics':ap.error('Different experiment directory; preserve the old run')
    if args.status:
        if not p:ap.error('No protocol')
        done,_=progress(p);print(f'WRITER {"ACTIVE" if b.writer_active(out) else "INACTIVE"}\nVERIFIED {done}/{len(p["scenarios"])*p["settings"]["reps"]} DATASETS');return 0
    if (args.resume or args.collect) and not p:ap.error('Cannot resume without protocol')
    if p and not(args.resume or args.collect):ap.error('Existing run: use --resume')
    if not p and out.exists() and any(f.name not in {'run.lock','preflight.log','environment.check.json'} for f in out.iterdir()):ap.error('Directory contains unrelated files')
    out.mkdir(parents=True,exist_ok=True)
    rs=shutil.which(args.rscript)
    if not rs:ap.error('Rscript not found')
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[k]='1'
    libs=b.libraries(args.r_library)
    if libs:env['PAPER1_CALIBRATION_LIBRARIES']=os.pathsep.join(libs)
    source={n:b.sha(ROOT/n) for n in FILES}
    try:
        with b.output_lock(out):
            print('NO DOWNLOADS; NO INSTALLS; ORIGINAL EXPERIMENTS UNCHANGED\nOUTPUT',out,flush=True)
            b.run([rs,'--vanilla',ROOT/'paired_count_calibration/worker.R','check',ROOT,out/'environment.check.json'],out/'preflight.log',env)
            observed=b.read(out/'environment.check.json');pin=b.read(ROOT/'inputs/R_versions.json')
            drift=observed['R']!=pin['R'] or any(observed['packages'][k]!=pin['packages'][k] for k in ('limma','edgeR'))
            if drift and not args.allow_version_drift:raise ValueError('R/limma/edgeR differ from pinned versions; see environment.check.json')
            if p:
                if p['output']!=str(out) or p['source']!=source:raise ValueError('Source or output changed')
                if any(p['environment'][k]!=observed[k] for k in ('R','packages','platform')):raise ValueError('Environment changed')
                for key in ('profile','reps','genes','seed'):
                    if getattr(args,key) is not None and getattr(args,key)!=p['settings'][key]:raise ValueError('Resume setting mismatch: '+key)
            else:
                profile=args.profile or 'full';defaults={'smoke':(1,250),'pilot':(10,500),'full':(500,3000)}
                reps=args.reps if args.reps is not None else defaults[profile][0]
                genes=args.genes if args.genes is not None else defaults[profile][1]
                seed=args.seed if args.seed is not None else 120261108
                if not 1<=reps<=99999 or genes<200 or not 100000000<=seed<=1000000000:raise ValueError('Require 1..99999 reps, >=200 genes, seed in 100000000..1000000000')
                p=dict(experiment='paired_tail_diagnostics',version=VERSION,output=str(out),source=source,
                    settings=dict(profile=profile,reps=reps,genes=genes,seed=seed),scenarios=scenarios(),
                    parameter_file=None,environment=observed,version_drift=drift,
                    controls=list(ARMS),selection='Six diagnostic scenarios selected after inspection of original grid; fresh replication seeds',
                    null='Global paired mean-contrast null; IID depth interpretation is marginal over resampled depths')
                p['run_id']=b.identity(p);b.atomic_json(out/'protocol.json',p);b.atomic_json(out/'environment.json',observed)
            if args.collect:report(p);return 0
            done,pending=progress(p);total=len(p['scenarios'])*p['settings']['reps']
            print(f'VERIFIED {done}/{total}; six arms, two tests; workers={args.workers}',flush=True)
            print(f'RESUME bash run_paired_tail_diagnostics.sh --resume --workers {args.workers} --output "{out}"'+(' --allow-version-drift' if drift else ''),flush=True)
            (out/'jobs').mkdir(exist_ok=True);jobs=[]
            for index,s,reps in pending:
                for j in range(0,len(reps),3):
                    path=out/'jobs'/f'{s["id"]}_{reps[j]:05d}.json';attempt=uuid.uuid4().hex
                    b.atomic_json(path,dict(protocol=str(out/'protocol.json'),scenario_index=index,replicates=reps[j:j+3],attempt=attempt))
                    jobs.append((path,s,reps[j:j+3],attempt))
            def execute(job):
                path,s,reps,attempt=job
                if any(b.sha(ROOT/n)!=h for n,h in source.items()):raise ValueError('Source changed while running')
                try:b.run([rs,'--vanilla',ROOT/'paired_tail_diagnostics/worker.R',ROOT,path],path.with_suffix('.log'),env)
                finally:
                    for rep in reps:
                        r=b.read(b.folder(p,s,rep)/'result.json')
                        if isinstance(r,dict) and r.get('attempt')==attempt:
                            try:receipt(p,s,rep)
                            except (OSError,ValueError,KeyError,TypeError):pass
                if not all(valid(p,s,r) for r in reps):raise ValueError('Incomplete results: '+str(path))
                return len(reps)
            pool=ThreadPoolExecutor(max_workers=args.workers);start=time.monotonic();new=0;last=0
            try:
                fs=[pool.submit(execute,j) for j in jobs]
                for f in as_completed(fs):
                    new+=f.result();now=time.monotonic()
                    if now-last>=10 or done+new==total:
                        print(f'TAIL DATASETS {done+new}/{total}; new {new}; ETA {(total-done-new)*(now-start)/new/3600:.2f} h',flush=True);last=now
            except BaseException:b.stop_children();raise
            finally:pool.shutdown(wait=True,cancel_futures=True)
            if any(b.sha(ROOT/n)!=h for n,h in source.items()):raise ValueError('Source changed; report withheld')
            report(p)
    except KeyboardInterrupt:
        b.stop_children();print('Interrupted; completed datasets retained. Resume the same output.');return 130
    except Exception as e:
        b.stop_children();print('STOP:',e,flush=True);return 1
    return 0
