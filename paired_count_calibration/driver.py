"""Independent datasets, reproducible seeds, validated receipts and compact reports."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time
import uuid
import zipfile
from . import VERSION

ROOT = Path(__file__).resolve().parent.parent
FILES = ('paper1_paired_count.py', 'run_paired_count_calibration.sh',
         'paired_count_calibration/__init__.py', 'paired_count_calibration/driver.py',
         'paired_count_calibration/engine.R', 'paired_count_calibration/worker.R',
         'paired_count_calibration/plot.R', 'analysis/common.R',
         'analysis/matched/functions.R', 'inputs/R_versions.json')
KEYS = {(w, t) for w in ('shared', 'design_specific') for t in ('model_based', 'HC3')}
CHILDREN, GUARD, STOP = set(), threading.Lock(), threading.Event()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read(path):
    try:
        return json.loads(Path(path).read_text(encoding='utf-8-sig'))
    except (OSError, ValueError):
        return None


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    os.replace(tmp, path)


def atomic_csv(path, rows):
    if not rows:
        raise ValueError('No report rows')
    tmp = path.with_suffix('.tmp')
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with tmp.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


@contextmanager
def output_lock(out):
    # OS locks disappear after a crash; never decide liveness from stale PID files.
    with (out / 'run.lock').open('a+b') as f:
        try:
            if os.name == 'posix':
                import fcntl
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            else:
                import msvcrt
                f.seek(0)
                if not f.read(1):
                    f.write(b'0'); f.flush()
                f.seek(0)
                msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as e:
            raise RuntimeError('Output has an active writer; do not delete run.lock') from e
        try:
            yield
        finally:
            if os.name == 'nt':
                f.seek(0); msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)


def writer_active(out):
    try:
        with output_lock(out):
            return False
    except RuntimeError:
        return True


def stop_children():
    STOP.set()
    with GUARD:
        children = list(CHILDREN)
    for p in children:
        if p.poll() is None:
            try:
                if os.name == 'posix':
                    os.killpg(p.pid, signal.SIGTERM)
                else:
                    p.terminate()
            except ProcessLookupError:
                pass
    for p in children:
        try:
            p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if os.name == 'posix':
                os.killpg(p.pid, signal.SIGKILL)
            else:
                p.kill()
            p.wait()


def run(cmd, log, env):
    with log.open('w', encoding='utf-8') as f:
        with GUARD:
            if STOP.is_set():
                raise RuntimeError('Interrupted')
            p = subprocess.Popen(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT,
                                 env=env, cwd=ROOT, start_new_session=os.name == 'posix')
            CHILDREN.add(p)
        try:
            code = p.wait()
            if code:
                raise RuntimeError(f'R failed ({code}): {log}\n' +
                                   log.read_text(encoding='utf-8', errors='replace')[-2500:])
        except BaseException:
            if p.poll() is None:
                p.terminate()
                try:
                    p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    p.kill(); p.wait()
            raise
        finally:
            with GUARD:
                CHILDREN.discard(p)


def scenarios(profile):
    result = []
    for i, (n, rho, depth, fraction) in enumerate(itertools.product(
            (40, 80, 160), (.3, .7), ('equal', 'iid_lognormal'), (0, .2)), 1):
        result.append(dict(id=f's{i:02d}', index=i, n=n, rho=rho, depth=depth,
                           exposure_fraction=fraction, design_seed=40000 + n))
    return [result[i] for i in (0, 3, 12, 15)] if profile == 'smoke' else result


def folder(p, s, rep):
    return Path(p['output']) / 'tasks' / s['id'] / f'r{rep:05d}'


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def validate_result(r, p, s, rep):
    if not isinstance(r, dict) or (r.get('run_id'), r.get('scenario'), r.get('replicate'), r.get('seed')) != (
            p['run_id'], s['id'], rep, p['settings']['seed'] + s['index'] * 100000 + rep):
        raise ValueError('Result identity mismatch')
    rows = r.get('metrics', [])
    if len(rows) != 4 or {(a.get('weight_mode'), a.get('test')) for a in rows} != KEYS:
        raise ValueError('Missing or duplicate test arm')
    for a in rows:
        g, k, n = a['genes'], a['retained'], a['available']
        if (any(type(v) is not int for v in (g, k, n)) or g != p['settings']['genes'] or
                not 0 <= n <= k <= g or k < 100 or a['filtered'] != g-k or
                a['unavailable_retained'] != k-n or a['paired_truth'] != 0):
            raise ValueError('Invalid availability totals')
        hist = a['p_histogram']
        if len(hist) != 20 or any(type(v) is not int or v < 0 for v in hist) or sum(hist) != n:
            raise ValueError('Invalid p-value histogram')
        if n == 0:
            raise ValueError('Entire arm unavailable; not a successful calibration dataset')
        for key in ('coverage95', 'delivered_coverage95', 'fpr_010', 'fpr_050', 'fpr_100'):
            if not finite(a[key]) or not 0 <= a[key] <= 1:
                raise ValueError('Invalid rate: '+key)
        for family in ('planned', 'retained'):
            for level in ('010', '050', '100'):
                count, fdp = a[f'discoveries_{family}_{level}'], a[f'fdp_{family}_{level}']
                if type(count) is not int or not 0 <= count <= n or fdp != int(count > 0):
                    raise ValueError('Invalid global-null FDP')
    for mode in ('shared', 'design_specific'):
        c = r['checks'][mode]
        if not finite(c['lmfit_max_difference']) or c['lmfit_max_difference'] >= 1e-7:
            raise ValueError('lmFit/kernel mismatch')
    return r


def valid_receipt(p, s, rep):
    d = folder(p, s, rep)
    receipt = read(d / 'receipt.json')
    if not isinstance(receipt, dict) or receipt.get('run_id') != p['run_id']:
        return False
    try:
        path = d / 'result.json'
        if receipt.get('sha256') != sha(path) or receipt.get('bytes') != path.stat().st_size:
            return False
        validate_result(read(path), p, s, rep)
        return True
    except (OSError, ValueError, TypeError, KeyError):
        return False


def publish_receipt(p, s, rep):
    d = folder(p, s, rep)
    path = d / 'result.json'
    validate_result(read(path), p, s, rep)
    atomic_json(d / 'receipt.json', dict(run_id=p['run_id'], sha256=sha(path), bytes=path.stat().st_size))


def progress(p):
    complete, pending = 0, []
    for i, s in enumerate(p['scenarios'], 1):
        missing = [r for r in range(1, p['settings']['reps']+1) if not valid_receipt(p, s, r)]
        complete += p['settings']['reps']-len(missing)
        pending.append((i, s, missing))
    return complete, pending


def mean_mc(values):
    values = [v for v in values if finite(v)]
    return (statistics.mean(values), statistics.stdev(values)/math.sqrt(len(values)) if len(values)>1 else None)


def wilson(events, total):
    z = 1.959963984540054
    v = events/total
    d = 1+z*z/total
    c = (v+z*z/(2*total))/d
    w = z*math.sqrt(v*(1-v)/total+z*z/(4*total*total))/d
    return (0 if events==0 else max(0, c-w)), (1 if events==total else min(1, c+w))


def collect(p, env, rscript):
    out = Path(p['output'])
    if any(sha(ROOT/n)!=digest for n,digest in p['source'].items()):
        raise ValueError('Source changed during execution; reports not published')
    complete, _ = progress(p)
    expected = len(p['scenarios'])*p['settings']['reps']
    if complete != expected:
        raise ValueError(f'Cannot report complete results: {complete}/{expected} verified')
    rows, warnings, histograms, grouped = [], [], [], {}
    for s in p['scenarios']:
        for rep in range(1, p['settings']['reps']+1):
            r = read(folder(p,s,rep)/'result.json')
            for message in r['warnings']:
                warnings.append(dict(scenario=s['id'],replicate=rep,message=message))
            for a in r['metrics']:
                row = {k:v for k,v in a.items() if k!='p_histogram'}
                row.update(scenario=s['id'],replicate=rep,seed=r['seed'],n=s['n'],rho=s['rho'],
                           depth=s['depth'],exposure_fraction=s['exposure_fraction'])
                rows.append(row)
                grouped.setdefault((s['id'],a['weight_mode'],a['test']),[]).append(a)
    summary = []
    for s in p['scenarios']:
        for mode, test in sorted(KEYS):
            a = grouped[(s['id'],mode,test)]
            row = {k:s[k] for k in ('n','rho','depth','exposure_fraction')}
            row.update(scenario=s['id'],weight_mode=mode,test=test,datasets=len(a),genes=p['settings']['genes'])
            for key in ('available','retained','coverage95','delivered_coverage95','mean_delta',
                        'rms_delta','mean_se','fpr_010','fpr_050','fpr_100'):
                row[key], row[key+'_mcse'] = mean_mc([x[key] for x in a])
            row['availability'] = row['available']/p['settings']['genes']
            for family in ('planned','retained'):
                for level in ('010','050','100'):
                    key = f'fdp_{family}_{level}'
                    events = sum(x[key] for x in a)
                    row[key] = events/len(a)
                    row[key+'_ci_low'], row[key+'_ci_high'] = wilson(events,len(a))
            summary.append(row)
            for j in range(20):
                histograms.append(dict(scenario=s['id'],weight_mode=mode,test=test,
                    lower=j*.05,upper=(j+1)*.05,count=sum(x['p_histogram'][j] for x in a)))
    atomic_csv(out/'summary.csv',summary)
    atomic_csv(out/'replicate_metrics.csv',rows)
    atomic_csv(out/'pvalue_histogram.csv',histograms)
    atomic_json(out/'warnings.json',warnings)
    run([rscript,'--vanilla',ROOT/'paired_count_calibration/plot.R',out],out/'plots.log',env)
    atomic_json(out/'completion.json',dict(run_id=p['run_id'],datasets=complete,
        status='complete',scope='global paired-contrast null; no non-null power or strong-FDR claim',
        validation_only=p['settings']['profile']!='full',warnings=len(warnings)))
    names = ['protocol.json','environment.json','completion.json','summary.csv','replicate_metrics.csv',
             'pvalue_histogram.csv','warnings.json','paired_FDR.pdf','paired_coverage.pdf','plots.log']
    if p['parameter_file']:
        names.append('parameters.csv')
    atomic_json(out/'review_manifest.json',{n:sha(out/n) for n in names})
    with zipfile.ZipFile(out/'review.zip.tmp','w',compression=zipfile.ZIP_DEFLATED) as z:
        for n in names+['review_manifest.json']:
            z.write(out/n,n)
        z.write(ROOT/'paired_count_calibration/README.md','METHODS.md')
        for n in FILES:
            z.write(ROOT/n,'source/'+n)
    os.replace(out/'review.zip.tmp',out/'review.zip')
    print('COMPLETE',complete,'datasets\nREVIEW',out/'review.zip',flush=True)


def libraries(explicit):
    found = list(explicit)
    for key in ('PAPER1_R_LIBRARIES','R_LIBS_USER','R_LIBS'):
        found.extend(p for p in os.environ.get(key,'').split(os.pathsep) if p)
    for base in (ROOT.parent,ROOT.parent/'paper1_public_evidence'):
        for pattern in ('paper1*/R_library','paper1*/runtime/R_library','p1_*/matched/R_library',
                        'p1_*/R_library','p1_*/runtime/R_library','p1_*/core/R_library',
                        'p1_*/matched/runtime/R_library','paper1_two_analyses/R_library'):
            found.extend(str(p) for p in base.glob(pattern))
    return list(dict.fromkeys(str(Path(p).resolve()) for p in found if Path(p).is_dir()))


def main(argv=None):
    def terminate(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, terminate)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--workers',type=int,default=36)
    ap.add_argument('--profile',choices=('smoke','pilot','full'))
    ap.add_argument('--reps',type=int)
    ap.add_argument('--genes',type=int)
    ap.add_argument('--seed',type=int)
    ap.add_argument('--parameters',type=Path,help='Optional mean_count,dispersion CSV; never old test results')
    ap.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
    ap.add_argument('--r-library',action='append',default=[])
    ap.add_argument('--allow-version-drift',action='store_true',help='Record non-pinned R/limma/edgeR; no silent fallback')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--status',action='store_true')
    ap.add_argument('--collect',action='store_true')
    args = ap.parse_args(argv)
    if args.workers<1:
        ap.error('--workers must be positive')
    if (args.resume or args.status or args.collect) and not args.output:
        ap.error('--output must be the actual existing results directory')
    out = (args.output or ROOT.parent/('p1_paircounts_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).expanduser().resolve()
    if out==ROOT or out.is_relative_to(ROOT) or ROOT.is_relative_to(out):
        ap.error('Output must be separate from source')
    p = read(out/'protocol.json')
    if args.status:
        if not p:
            ap.error('No protocol.json in output')
        n,_=progress(p);total=len(p['scenarios'])*p['settings']['reps']
        print(f'WRITER {"ACTIVE" if writer_active(out) else "INACTIVE"}\nVERIFIED DATASETS {n}/{total} ({100*n/total:.1f}%)\nREMAINING {total-n}')
        return 0
    if not p and (args.resume or args.collect):
        ap.error('Cannot resume/collect without protocol.json')
    if p and not (args.resume or args.collect):
        ap.error('Existing protocol: use --resume or a new output')
    allowed_preflight={'run.lock','preflight.log','environment.check.json'}
    if not p and out.exists() and any(f.name not in allowed_preflight for f in out.iterdir()):
        ap.error('Output contains unrelated files; choose a new directory')
    out.mkdir(parents=True,exist_ok=True)
    rscript = shutil.which(args.rscript)
    if not rscript:
        ap.error('Rscript not found; activate the existing R environment or use --rscript')
    env = os.environ.copy()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[k]='1'
    libs=libraries(args.r_library)
    if libs:
        env['PAPER1_CALIBRATION_LIBRARIES']=os.pathsep.join(libs)
    source={n:sha(ROOT/n) for n in FILES}
    try:
        with output_lock(out):
            print('NO DOWNLOADS; NO PACKAGE INSTALLS; ORIGINAL RESULTS UNCHANGED\nOUTPUT',out,flush=True)
            run([rscript,'--vanilla',ROOT/'paired_count_calibration/worker.R','check',ROOT,out/'environment.check.json'],
                out/'preflight.log',env)
            observed=read(out/'environment.check.json')
            pin=read(ROOT/'inputs/R_versions.json')
            mismatch=(observed['R']!=pin['R'] or any(observed['packages'][k]!=pin['packages'][k] for k in ('limma','edgeR')))
            if mismatch and not args.allow_version_drift:
                raise ValueError('R/limma/edgeR differ from manuscript pins. See environment.check.json; use the pinned environment or explicitly --allow-version-drift.')
            if p:
                if p['output']!=str(out) or p['source']!=source:
                    raise ValueError('Source or output changed; do not mix incompatible experiments')
                for key in ('profile','reps','genes','seed'):
                    if getattr(args,key) is not None and getattr(args,key)!=p['settings'][key]:
                        raise ValueError('Resume setting mismatch: '+key)
                old=read(out/'environment.json')
                if any(old[k]!=observed[k] for k in ('R','packages','platform')):
                    raise ValueError('R environment differs from the original run')
                if p['parameter_file'] and sha(p['parameter_file'])!=p['parameter_sha256']:
                    raise ValueError('Parameter input changed')
                if args.parameters and sha(args.parameters)!=p['parameter_sha256']:
                    raise ValueError('Resume parameter file mismatch')
            else:
                profile=args.profile or 'full'
                defaults={'smoke':(1,250),'pilot':(20,1000),'full':(1000,3000)}
                reps=args.reps if args.reps is not None else defaults[profile][0]
                genes=args.genes if args.genes is not None else defaults[profile][1]
                seed=args.seed if args.seed is not None else 20261107
                if not 1<=reps<=99999 or genes<200 or not 1<=seed<=2_000_000_000:
                    raise ValueError('Need 1..99999 reps, >=200 genes, seed 1..2000000000')
                param_file=None
                if args.parameters:
                    if not args.parameters.is_file():
                        raise ValueError('Parameter file not found')
                    param_file=str(out/'parameters.csv')
                    shutil.copy2(args.parameters,param_file)
                p=dict(version=VERSION,output=str(out),source=source,
                       settings=dict(profile=profile,reps=reps,genes=genes,seed=seed),
                       scenarios=scenarios(profile),parameter_file=param_file,
                       parameter_sha256=sha(param_file) if param_file else None,
                       normalization='filterByExpr(X0) -> TMM -> voom adaptive.span=TRUE',
                       truth='all paired mean contrasts are zero by within-exposure exchangeability',
                       iid_depth_scope='marginal over independently resampled library depth',
                       full_plan_default=(profile=='full' and reps==1000 and genes==3000),
                       environment=observed,version_drift=mismatch)
                p['run_id']=identity(p)
                atomic_json(out/'environment.json',observed)
                atomic_json(out/'protocol.json',p)
            if args.collect:
                collect(p,env,rscript)
                return 0
            done,pending=progress(p)
            total=len(p['scenarios'])*p['settings']['reps']
            print(f'VERIFIED {done}/{total}; workers={args.workers}; genes={p["settings"]["genes"]}',flush=True)
            print(f'RESUME bash run_paired_count_calibration.sh --resume --output "{out}" --workers {args.workers}'+
                  (' --allow-version-drift' if mismatch else ''),flush=True)
            batches=[]
            (out/'jobs').mkdir(exist_ok=True)
            for index,s,reps in pending:
                for j in range(0,len(reps),5):
                    batch=dict(protocol=str(out/'protocol.json'),scenario_index=index,
                               replicates=reps[j:j+5],attempt=uuid.uuid4().hex)
                    name=f'{s["id"]}_{reps[j]:05d}'
                    path=out/'jobs'/f'{name}.json'
                    atomic_json(path,batch)
                    batches.append((path,s,batch['replicates'],batch['attempt']))
            def execute(item):
                path,s,reps,attempt=item
                if any(sha(ROOT/n)!=digest for n,digest in p['source'].items()):
                    raise ValueError('Source changed during execution; refusing a mixed-source batch')
                try:
                    run([rscript,'--vanilla',ROOT/'paired_count_calibration/worker.R','run',ROOT,path],path.with_suffix('.log'),env)
                finally:
                    # Completed datasets survive a later failure in the same five-replicate batch.
                    for rep in reps:
                        candidate=read(folder(p,s,rep)/'result.json')
                        if isinstance(candidate,dict) and candidate.get('attempt')==attempt:
                            try:
                                publish_receipt(p,s,rep)
                            except (OSError,ValueError,KeyError,TypeError):
                                pass
                if not all(valid_receipt(p,s,r) for r in reps):
                    raise ValueError('Incomplete dataset in '+str(path))
                return len(reps)
            start=time.monotonic();new=0;last_print=0
            pool=ThreadPoolExecutor(max_workers=args.workers)
            try:
                futures=[pool.submit(execute,b) for b in batches]
                for future in as_completed(futures):
                    new+=future.result()
                    elapsed=time.monotonic()-start
                    eta=(total-done-new)*elapsed/new/3600
                    if time.monotonic()-last_print>=10 or done+new==total:
                        print(f'PAIRED COUNT DATASETS {done+new}/{total}; new {new}; ETA {eta:.2f} h',flush=True)
                        last_print=time.monotonic()
            except BaseException:
                stop_children()
                raise
            finally:
                pool.shutdown(wait=True,cancel_futures=True)
            collect(p,env,rscript)
    except KeyboardInterrupt:
        stop_children()
        print('Interrupted. Completed datasets are retained; resume the printed output directory.',file=sys.stderr)
        return 130
    except Exception as e:
        stop_children()
        print('STOP:',e,file=sys.stderr)
        return 1
    return 0
