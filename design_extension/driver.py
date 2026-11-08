"""Resumable extension, separate from original analyses and publication outputs."""
import argparse
import concurrent.futures as cf
import contextlib
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
import zipfile
import numpy as np
import scipy
from . import VERSION
from .study import create_designs, scenarios, read_design, write_csv, MODELS
from .geometry import design_diagnostics
from .report import build_report, csv_rows, COUNT_MODELS, number

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')
CHILDREN, GUARD, STOP = set(), threading.Lock(), threading.Event()
COUNT_FILES = ('metrics.csv','failure_counts.csv','diagnostics.csv',
               'retry_attempts.csv','messages.csv','unavailable_genes.csv')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def identity(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def read_json(path):
    try:
        value=json.loads(Path(path).read_text(encoding='utf-8-sig'))
        return value if isinstance(value,dict) else {}
    except (OSError,ValueError):return {}


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
    os.replace(tmp,path)


@contextlib.contextmanager
def lock_output(out):
    with (out/'run.lock').open('a+b') as f:
        try:
            if os.name=='posix':
                import fcntl
                fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
            else:
                import msvcrt
                f.seek(0)
                if not f.read(1):f.write(b'0');f.flush()
                f.seek(0);msvcrt.locking(f.fileno(),msvcrt.LK_NBLCK,1)
        except OSError as e:raise RuntimeError('Output has an active writer. Do not delete run.lock.') from e
        try:yield
        finally:
            if os.name=='nt':f.seek(0);msvcrt.locking(f.fileno(),msvcrt.LK_UNLCK,1)


def writer_active(out):
    if not (out/'run.lock').exists():return False
    try:
        with lock_output(out):return False
    except RuntimeError:return True


def stop_children():
    STOP.set()
    with GUARD:children=list(CHILDREN)
    for p in children:
        if p.poll() is None:
            try:
                if os.name=='posix':os.killpg(p.pid,signal.SIGTERM)
                else:p.terminate()
            except ProcessLookupError:pass
    for p in children:
        try:p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                if os.name=='posix':os.killpg(p.pid,signal.SIGKILL)
                else:p.kill()
            except ProcessLookupError:pass
            p.wait()


def command(cmd, env, log, cwd):
    log.parent.mkdir(parents=True,exist_ok=True)
    with GUARD:
        if STOP.is_set():raise RuntimeError('Interrupted')
        f=log.open('w',encoding='utf-8')
        p=subprocess.Popen(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,env=env,cwd=cwd,
                           start_new_session=os.name=='posix')
        CHILDREN.add(p)
    try:
        if p.wait():raise RuntimeError('Task failed; see '+str(log)+'\n'+log.read_text(encoding='utf-8',errors='replace')[-2200:])
    except BaseException:
        stop_children()
        raise
    finally:
        f.close()
        with GUARD:CHILDREN.discard(p)


INPUTS=('prepared/mouse_units.csv','prepared/mouse_whole_counts.csv.gz','prepared/mouse_cortex_counts.csv.gz')


def discover_core(explicit, parent):
    if explicit:
        p=Path(explicit).expanduser().resolve()
        if all((p/n).is_file() for n in INPUTS):return p
        if all((p/'core'/n).is_file() for n in INPUTS):return p/'core'
        raise ValueError('Need prepared mouse counts AND metadata in --data-work; review-only exports are insufficient')
    found=[]
    for pattern in ('p1_run_*/core','paper1_complete_from_public*/core','paper1_from_public*/core','paper1_from_raw*/core'):
        for p in parent.glob(pattern):
            if all((p/n).is_file() for n in INPUTS):found.append(p)
    if not found:raise ValueError('Prepared public counts not found. Use --data-work /path/to/completed/core')
    return max(found,key=lambda p:(p/'prepared/mouse_units.csv').stat().st_mtime).resolve()


def libraries(parent, explicit):
    paths=[Path(p).expanduser() for p in explicit]
    for key in ('PAPER1_R_LIBRARIES','R_LIBS','R_LIBS_USER'):
        paths += [Path(p) for p in os.environ.get(key,'').split(os.pathsep) if p]
    for base in (parent,parent/'paper1_public_evidence'):
        if base.is_dir():
            for child in sorted(base.iterdir()):
                if child.is_dir() and child.name.startswith(('paper1','p1_')):
                    paths += [child/'R_library',child/'runtime/R_library',child/'matched/runtime/R_library']
    paths += [parent/'paper1_public_evidence/paper1_two_analyses/R_library']
    return list(dict.fromkeys(str(p.resolve()) for p in paths if p.is_dir()))


def receipt_problem(out,job,run_id):
    folder=out/'tasks'/job['id'];r=read_json(folder/'receipt.json')
    if not (folder/'receipt.json').is_file():return 'missing receipt'
    if not r:return 'unreadable or non-object receipt'
    if r.get('run_id')!=run_id or r.get('job_hash')!=identity(job):return 'receipt identity mismatch'
    if not isinstance(r.get('files'),dict):return 'invalid file manifest'
    required=set(COUNT_FILES) if job['phase']=='counts' else {'metrics.csv','sensitivity.csv'}
    if set(r['files'])!=required:return 'incomplete file manifest'
    for name,record in r['files'].items():
        if not isinstance(record,dict) or set(record)!={'bytes','sha256'}:return 'invalid file record'
        if type(record['bytes']) is not int or record['bytes']<0 or not isinstance(record['sha256'],str):
            return 'invalid size or hash type'
        path=(folder/name).resolve()
        try:
            if not path.is_relative_to(folder.resolve()) or not path.is_file():return 'missing or unsafe file: '+name
            if path.stat().st_size!=record['bytes'] or sha(path)!=record['sha256']:return 'file checksum mismatch: '+name
        except OSError:return 'unreadable result file: '+name
    return None


def valid_receipt(out,job,run_id):
    return receipt_problem(out,job,run_id) is None


def execution_identity(protocol):
    return identity({k:protocol[k] for k in ('output','data_work','designs','scenarios','jobs','design_hashes','empirical_design')})


def check_output_location(out, code_root, data_work=None):
    out=Path(out).resolve();code_root=Path(code_root).resolve()
    for protected in (code_root, Path(data_work).resolve() if data_work else None):
        if protected and (out==protected or out.is_relative_to(protected) or protected.is_relative_to(out)):
            raise ValueError('Output must not overlap code or prepared inputs')
    allowed={'run.lock','R_check.log','environment.check.json'}
    if out.exists() and not (out/'protocol.json').exists() and any(p.name not in allowed for p in out.iterdir()):
        raise ValueError('Output contains unrelated files; choose a new directory')


def validate_metrics(rows, job):
    models = COUNT_MODELS if job['phase']=='counts' else MODELS
    expected = {(m, r) for m in models for r in range(job['start'],job['end'])}
    keys = {(r['model'],int(r['replicate'])) for r in rows}
    if len(rows)!=len(expected) or keys!=expected:
        raise ValueError('Incomplete, duplicate or unexpected output in '+job['id'])
    for row in rows:
        if (row['design']!=job['design'] or row['split']!=job['split'] or
            int(row['scenario'])!=job['scenario']['scenario']):
            raise ValueError('Output identity mismatch in '+job['id'])
        availability = number(row['availability'] if job['phase']=='counts' else row['available'])
        if availability is None or not 0<=availability<=1:
            raise ValueError('Invalid availability in '+job['id'])
        for metric in ('FDP','coverage','power'):
            v=number(row.get(metric))
            if v is not None and not 0<=v<=1:
                raise ValueError('Invalid '+metric+' in '+job['id'])
        if availability>0 and (number(row.get('FDP')) is None or number(row.get('coverage')) is None):
            raise ValueError('Missing performance for available inference in '+job['id'])
        if availability==0 and number(row.get('FDP')) is not None:
            raise ValueError('Unavailable inference reported as a usable FDP in '+job['id'])
        if job['phase']=='counts':
            g=int(row['genes']); n=int(row['available']); k=int(row['filtered'])
            discoveries=int(row['discoveries']); false=int(row['false_discoveries']); true=int(row['true_discoveries'])
            if (g<=0 or not 0<=n<=g-k or not 0<=k<=g or
                abs(n/g-availability)>1e-10 or discoveries!=false+true or
                not 0<=false<=int(row['null_genes']) or not 0<=true<=int(row['alternative_genes']) or
                int(row['null_genes'])+int(row['alternative_genes'])!=g):
                raise ValueError('Inconsistent count-performance totals in '+job['id'])


def validate_count_health(rows, job):
    _,x,g,c=read_design(job['design_path'])
    geometry=design_diagnostics(x,g,c)
    required=['baseline','cage_correlation','dream_native']
    if geometry['within_estimable'] and geometry['fixed_residual_df']>=5:
        required+=['cage_fixed','mundlak_native']
    failed=[m for m in required if all(float(r['availability'])==0 for r in rows if r['model']==m)]
    if failed:
        raise ValueError('No usable inference in an entire batch for '+', '.join(failed)+
                         '; inspect '+job['output']+'/run.log before resuming')


def execute(job,protocol,env):
    out=Path(protocol['output']);folder=out/'tasks'/job['id'];folder.mkdir(parents=True,exist_ok=True)
    atomic_json(folder/'job.json',job)
    started=time.time();source=out/'source'
    if job['phase']=='counts':
        cmd=[protocol['rscript'],'--vanilla',source/'design_extension/count_bridge.R','job',folder/'job.json']
    else:cmd=[sys.executable,'-B','-m','design_extension.worker',folder/'job.json']
    task_env=env.copy()
    task_env['PAPER1_RUN_ID']=protocol['run_id']
    task_env['PAPER1_JOB_HASH']=identity(job)
    command(cmd,task_env,folder/'run.log',source)
    rows=csv_rows(folder/'metrics.csv')
    validate_metrics(rows,job)
    if job['phase']=='counts':validate_count_health(rows,job)
    files=COUNT_FILES if job['phase']=='counts' else ['metrics.csv','sensitivity.csv']
    receipt=dict(run_id=protocol['run_id'],job_hash=identity(job),seconds=time.time()-started,
                 work_items=job['end']-job['start'],files={n:dict(bytes=(folder/n).stat().st_size,sha256=sha(folder/n)) for n in files})
    atomic_json(folder/'receipt.json',receipt)
    return receipt


def current_partial(out,jobs,done_ids,protocol,state):
    count=0
    for j in jobs:
        if j['id'] in done_ids:continue
        p=read_json(out/'tasks'/j['id']/'progress.json')
        if (not state.get('attempt_id') or p.get('attempt_id')!=state['attempt_id'] or
            p.get('run_id')!=protocol['run_id'] or p.get('job_hash')!=identity(j)):
            continue
        done=p.get('done')
        if type(done) is int and 0<=done<=j['end']-j['start']:count+=done
    return count


def status(out):
    p=read_json(out/'protocol.json')
    if not p:print('No protocol. Check environment/input preparation log.');return
    if p.get('version',VERSION)!=VERSION:
        raise ValueError('Use the original code for legacy run status; this version has different checkpoint files')
    active=writer_active(out);state=read_json(out/'state.json')
    print('STATE',state.get('status'),'WRITER', 'ACTIVE' if active else 'INACTIVE')
    for phase in ('gaussian','counts'):
        jobs=[j for j in p['jobs'] if j['phase']==phase]
        done=[j for j in jobs if valid_receipt(out,j,p['run_id'])]
        total=sum(j['end']-j['start'] for j in jobs);completed=sum(j['end']-j['start'] for j in done)
        print(phase,f'{len(done)}/{len(jobs)} batches; {completed}/{total} datasets')
        if active and phase=='counts':
            finished={j['id'] for j in done}
            partial=current_partial(out,jobs,finished,p,state)
            print(f'  Reported datasets in unfinished batches: {partial} (not yet receipt-verified).')
        if active and len(done)>=3 and completed and completed<total:
            seconds=sum(read_json(out/'tasks'/j['id']/'receipt.json')['seconds'] for j in done)
            hours=seconds/completed*(total-completed)/max(1,state.get('workers',1))/3600
            print(f'  Rough phase-only remaining time: {hours:.1f} h; pilot/full sizes and contention matter.')
    if not active:print('No active writer. Resume requires --resume and the same output/profile/settings.')


def make_jobs(designs,scens,out,seed,phase,gaussian_reps,count_reps,count_genes):
    jobs=[]
    for di,d in enumerate(designs):
        for s in scens:
            # A count dataset is the smallest resumable statistical unit.
            for kind,reps,batch,genes in [('gaussian',gaussian_reps,100,256),('counts',count_reps,1,count_genes)]:
                if kind=='counts' and phase=='gaussian':continue
                for start in range(0,reps,batch):
                    identifier=f'{kind}_{d["design"]}_s{s["scenario"]}_{start:05d}'
                    jobs.append(dict(id=identifier,phase=kind,design=d['design'],split=d['split'],design_index=di,
                        design_path=d['path'],scenario=s,start=start,end=min(start+batch,reps),genes=genes,seed=seed,
                        output=str(out/'tasks'/identifier),config=str(out/'config.json')))
    return jobs


def collect(out):
    p=read_json(out/'protocol.json')
    if writer_active(out):raise RuntimeError('Wait for the writer before collecting the review package')
    if not p or not (out/'review/REPORT.md').exists():raise ValueError('No finished review report')
    if not all(valid_receipt(out,j,p['run_id']) for j in p['jobs']):raise ValueError('Incomplete or damaged run')
    target=out.parent/(out.name+'_review.zip')
    if target.exists():target=out.parent/(out.name+'_review_'+datetime.now().strftime('%H%M%S')+'.zip')
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        files=[out/'protocol.json',out/'state.json',out/'environment.json',out/'calibration_info.json']
        files+=list((out/'review').glob('*'))
        for path in files:
            if path.is_file():z.write(path,arcname=out.name+'/'+str(path.relative_to(out)).replace('\\','/'))
    print('REVIEW PACKAGE',target)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    p.add_argument('--data-work',type=Path)
    p.add_argument('--workers',type=int,default=36)
    p.add_argument('--profile',choices=('pilot','full'),default='pilot')
    p.add_argument('--phase',choices=('all','gaussian'),default='all')
    p.add_argument('--count-reps',type=int)
    p.add_argument('--count-genes',type=int)
    p.add_argument('--gaussian-reps',type=int)
    p.add_argument('--seed',type=int,default=20261101)
    p.add_argument('--fdr-screen',type=float,default=.075)
    p.add_argument('--rscript')
    p.add_argument('--r-library',action='append',default=[])
    p.add_argument('--resume',action='store_true')
    p.add_argument('--status',action='store_true')
    p.add_argument('--collect',action='store_true')
    p.add_argument('--design-csv',type=Path,help='Inspect metadata only; requires animal_id/genotype/sex/age/cage columns')
    return p


def main(argv=None):
    def interrupted(signum, frame):
        raise KeyboardInterrupt('Received signal '+str(signum))
    signal.signal(signal.SIGTERM,interrupted)
    STOP.clear()
    a=parser().parse_args(argv)
    if a.design_csv:
        _,x,g,c=read_design(a.design_csv)
        print(json.dumps(design_diagnostics(x,g,c),indent=2,allow_nan=False));return 0
    if (a.resume or a.status or a.collect) and a.output is None:raise ValueError('Specify the exact --output')
    if a.status:status(a.output.resolve());return 0
    if a.collect:collect(a.output.resolve());return 0
    cpu=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1
    if not 1<=a.workers<=cpu:raise ValueError(f'workers must be in 1..{cpu}')
    if a.seed<0 or a.seed>100000000 or not 0<a.fdr_screen<1:raise ValueError('Invalid seed/risk screen')
    parent=a.output.expanduser().resolve().parent if a.output is not None else HERE.parent.parent
    out=(a.output or parent/('p1_design_fixed_'+datetime.now().strftime('%m%d_%H%M%S'))).resolve()
    if (out/'protocol.json').is_file() and read_json(out/'protocol.json').get('version')!=VERSION:
        raise ValueError('This source cannot resume a legacy run. Use the original source to resume it, or repair its exported report without fitting')
    selected_input=discover_core(a.data_work,parent) if a.phase=='all' else None
    check_output_location(out,HERE.parent,selected_input)
    out.mkdir(parents=True,exist_ok=True)
    with lock_output(out):
        old=read_json(out/'protocol.json')
        if old and not a.resume:raise ValueError('Run exists; use --resume with the same settings')
        count_reps=a.count_reps if a.count_reps is not None else (3 if a.profile=='pilot' else 50)
        count_genes=a.count_genes if a.count_genes is not None else (600 if a.profile=='pilot' else 0)
        gaussian_reps=a.gaussian_reps if a.gaussian_reps is not None else (100 if a.profile=='pilot' else 500)
        if min(count_reps,gaussian_reps)<1 or count_genes<0 or 0<count_genes<100:raise ValueError('Invalid simulation sizes')
        if a.seed+11*10007+100003*(count_reps-1)>=2147483647:
            raise ValueError('Count repetition/seed range exceeds the R generator integer limit')
        env=os.environ.copy();env.update({k:'1' for k in THREADS})
        env['BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS']='FALSE'
        code_files=sorted(p for p in HERE.rglob('*') if p.suffix in ('.py','.R') and '__pycache__' not in p.parts)
        code={str(p.relative_to(HERE)):sha(p) for p in code_files}
        data_work=None;inputs={};r_info={};rscript=None;empirical_design=None
        if a.phase=='all':
            data_work=selected_input
            for n in INPUTS:print('HASH INPUT',data_work/n,flush=True);inputs[n]=sha(data_work/n)
            _,x,groups,covariates=read_design(data_work/'prepared/mouse_units.csv')
            empirical_design=dict(source='prepared mouse metadata',**design_diagnostics(x,groups,covariates))
            receipt=read_json(data_work/'stages/prepare_mouse.json')
            if receipt:
                for n in INPUTS:
                    if receipt.get('outputs',{}).get(n)!=inputs[n]:raise ValueError('Prepared-input receipt mismatch: '+n)
            rscript=a.rscript or shutil.which('Rscript')
            if not rscript:raise ValueError('Activate paper1 or pass --rscript; no packages are installed')
            if os.name=='nt' and (Path(rscript).parent/'x64/Rscript.exe').is_file():rscript=str(Path(rscript).parent/'x64/Rscript.exe')
            env['PAPER1_R_LIBRARIES']=os.pathsep.join(libraries(parent,a.r_library))
            command([rscript,'--vanilla',HERE/'count_bridge.R','check',out/'environment.check.json'],env,out/'R_check.log',HERE.parent)
            r_info=read_json(out/'environment.check.json')
        runtime=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,R=r_info)
        scientific=dict(version=VERSION,profile=a.profile,phase=a.phase,count_reps=count_reps,count_genes=count_genes,
                        gaussian_reps=gaussian_reps,seed=a.seed,fdr_screen=a.fdr_screen,code=code,inputs=inputs,runtime=runtime)
        run_id=identity(scientific)
        if old and old['run_id']!=run_id:raise ValueError('Source/input/environment/settings changed; use a new output directory')
        if old and (old['output']!=str(out) or old['data_work']!=(str(data_work) if data_work else None)):
            raise ValueError('Run or input location changed; use original locations or a fresh run')
        atomic_json(out/'environment.json',runtime)
        if not old:
            for path in code_files:
                dest=out/'source/design_extension'/path.relative_to(HERE);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,dest)
            designs=create_designs(out/'designs',a.seed);scens=scenarios(a.profile)
            jobs=make_jobs(designs,scens,out,a.seed,a.phase,gaussian_reps,count_reps,count_genes)
            p=dict(**scientific,run_id=run_id,output=str(out),data_work=str(data_work) if data_work else None,
                   rscript=rscript,designs=designs,scenarios=scens,jobs=jobs,
                   empirical_design=empirical_design,
                   design_hashes={d['path']:sha(d['path']) for d in designs})
            p['execution_manifest_hash']=execution_identity(p)
            atomic_json(out/'protocol.json',p)
            atomic_json(out/'config.json',dict(output=str(out),data_work=str(data_work),count_genes=count_genes,seed=a.seed))
        else:p=old
        if p.get('execution_manifest_hash')!=execution_identity(p):
            raise ValueError('Execution manifest changed')
        expected_config=dict(output=str(out),data_work=str(data_work),count_genes=count_genes,seed=a.seed)
        if read_json(out/'config.json')!=expected_config:
            raise ValueError('Run configuration differs from its signed protocol')
        for path,digest in p['design_hashes'].items():
            if sha(path)!=digest:raise ValueError('Design metadata changed: '+path)
        for name,digest in code.items():
            if sha(out/'source/design_extension'/name)!=digest:raise ValueError('Frozen source changed: '+name)
        if a.phase=='all':
            prep=read_json(out/'calibration_receipt.json')
            names=('calibration.rds','calibration_parameters.csv','calibration_info.json')
            valid=prep.get('run_id')==run_id and all((out/n).is_file() and sha(out/n)==prep.get('files',{}).get(n) for n in names)
            if not valid:
                print('FRESH COUNT CALIBRATION',flush=True)
                command([rscript,'--vanilla',out/'source/design_extension/count_bridge.R','prepare',out/'config.json'],env,out/'prepare.log',out/'source')
                atomic_json(out/'calibration_receipt.json',dict(run_id=run_id,files={n:sha(out/n) for n in names}))
        pending=[j for j in p['jobs'] if not valid_receipt(out,j,run_id)]
        damaged=[dict(task=j['id'],reason=receipt_problem(out,j,run_id)) for j in pending
                 if (out/'tasks'/j['id']/'receipt.json').exists()]
        if damaged:
            atomic_json(out/'receipt_audit.json',dict(invalid_receipts=damaged))
            print('INVALID RECEIPTS',len(damaged),'queued for recomputation; see receipt_audit.json',flush=True)
        print('OUTPUT',out,'WORKERS',a.workers,'PENDING BATCHES',len(pending),flush=True)
        attempt_id=uuid.uuid4().hex
        env['PAPER1_ATTEMPT_ID']=attempt_id
        atomic_json(out/'state.json',dict(status='RUNNING',pid=os.getpid(),workers=a.workers,
                                       started=time.time(),attempt_id=attempt_id))
        try:
            for phase in ('gaussian','counts'):
                phase_jobs=[j for j in pending if j['phase']==phase]
                pool=cf.ThreadPoolExecutor(max_workers=a.workers)
                try:
                    futures={pool.submit(execute,j,p,env):j for j in phase_jobs}
                    for i,f in enumerate(cf.as_completed(futures),1):
                        result=f.result()
                        print(phase,f'{i}/{len(phase_jobs)} new batches',f'{result["seconds"]:.1f}s',flush=True)
                except BaseException:
                    stop_children()
                    for f in futures:f.cancel()
                    raise
                finally:pool.shutdown(wait=True,cancel_futures=True)
            if not all(valid_receipt(out,j,run_id) for j in p['jobs']):raise ValueError('Final receipt verification failed')
            (out/'review').mkdir(exist_ok=True)
            issues=build_report(out,p)
            atomic_json(out/'state.json',dict(status='COMPLETED_WITH_FAILURES' if issues else 'COMPLETED',
                workers=a.workers,ended=time.time(),whole_arm_failures=issues,profile=a.profile))
        except BaseException as e:
            stop_children();atomic_json(out/'state.json',dict(status='INTERRUPTED' if isinstance(e,KeyboardInterrupt) else 'FAILED',
                workers=a.workers,error=str(e),ended=time.time()));raise
    print('REPORT',out/'review/report.html',flush=True)
    collect(out)
    return 0
