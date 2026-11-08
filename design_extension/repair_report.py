"""Rebuild reports from completed exported metrics without fitting any models."""
import argparse
from collections import Counter
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile

import numpy as np
import scipy

from . import VERSION
from .driver import identity, read_json, sha, validate_metrics, atomic_json
from .report import build_report, csv_rows, COUNT_MODELS, number
from .study import MODELS, write_csv


SCIENCE = ('version','profile','phase','count_reps','count_genes','gaussian_reps',
           'seed','fdr_screen','code','inputs','runtime')
EXECUTION = ('output','data_work','designs','scenarios','jobs','design_hashes','empirical_design')


@contextmanager
def input_read_lock(source):
    lock=source/'run.lock'
    if lock.is_file() and os.name=='posix':
        import fcntl
        with lock.open('rb') as stream:
            try:fcntl.flock(stream,fcntl.LOCK_SH|fcntl.LOCK_NB)
            except BlockingIOError as e:raise ValueError('Source has an active writer; wait for completion') from e
            try:yield
            finally:fcntl.flock(stream,fcntl.LOCK_UN)
    else:yield


def load_source(source):
    p=read_json(source/'protocol.json')
    if not set(SCIENCE+EXECUTION).issubset(p):raise ValueError('Not a complete design-extension protocol')
    if p.get('run_id')!=identity({k:p[k] for k in SCIENCE}):raise ValueError('Source scientific signature mismatch')
    if p.get('execution_manifest_hash')!=identity({k:p[k] for k in EXECUTION}):raise ValueError('Source execution signature mismatch')
    if p['version'] not in ('0.1.0',VERSION):raise ValueError('Unsupported source protocol version')
    if read_json(source/'state.json').get('status') not in ('COMPLETED','COMPLETED_WITH_FAILURES'):
        raise ValueError('Report repair requires a completed run or its completed review export')
    if (source/'environment.json').is_file() and read_json(source/'environment.json')!=p['runtime']:
        raise ValueError('Source environment file does not match the protocol')
    if p['profile'] not in ('pilot','full') or p['phase'] not in ('all','gaussian'):
        raise ValueError('Invalid source profile or phase')
    if not {'driver.py','count_bridge.R'}.issubset(p['code']):raise ValueError('Wrong study source manifest')
    if len({j['id'] for j in p['jobs']})!=len(p['jobs']):raise ValueError('Duplicate source task IDs')
    return p


def snapshot(source,paths):
    return {str(p.relative_to(source)).replace('\\','/'):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in paths}


def validate_export(rows,jobs,phase,protocol,genes=None):
    models=COUNT_MODELS if phase=='counts' else MODELS
    owners={}
    for j in jobs:
        if j['phase']!=phase:continue
        if not 0<=j['start']<j['end']:raise ValueError('Invalid replicate interval')
        for r in range(j['start'],j['end']):
            for m in models:
                key=(j['design'],int(j['scenario']['scenario']),r,m)
                if key in owners:raise ValueError('Overlapping source jobs')
                owners[key]=j['id']
    grouped={j['id']:[] for j in jobs if j['phase']==phase}
    seen=set()
    scenarios={int(s['scenario']):s for s in protocol['scenarios']}
    for row in rows:
        key=(row['design'],int(row['scenario']),int(row['replicate']),row['model'])
        if key in seen:raise ValueError('Duplicate dataset/model record: '+str(key))
        if key not in owners:raise ValueError('Unexpected dataset/model record: '+str(key))
        seen.add(key);grouped[owners[key]].append(row)
        if phase=='counts':
            for k in ('cage_log2_sd','cage_mean_amplitude','de_fraction'):
                if number(row.get(k))!=float(scenarios[key[1]][k]):raise ValueError('Scenario parameter mismatch: '+k)
            if genes is not None and int(row['genes'])!=genes:raise ValueError('Count gene background mismatch')
            n=int(row['available']);total=int(row['genes']);alt=int(row['alternative_genes'])
            false=int(row['false_discoveries']);true=int(row['true_discoveries']);disc=int(row['discoveries'])
            if disc>n:raise ValueError('More discoveries than available tests')
            expected=false/max(1,disc) if n else None
            actual=number(row.get('FDP'))
            if (actual is None)!=(expected is None) or (actual is not None and not math.isclose(actual,expected,abs_tol=1e-12)):
                raise ValueError('FDP does not match discovery totals')
            power=true/alt if n and alt else None
            actual=number(row.get('power'))
            if (actual is None)!=(power is None) or (actual is not None and not math.isclose(actual,power,abs_tol=1e-12)):
                raise ValueError('Power does not match discovery totals')
            if not math.isclose(float(row['availability']),n/total,abs_tol=1e-12):raise ValueError('Availability mismatch')
    if seen!=set(owners):raise ValueError(f'{phase}: missing {len(set(owners)-seen)} dataset/model records')
    for j in jobs:
        if j['phase']==phase:validate_metrics(grouped[j['id']],j)
    return dict(planned_rows=len(owners),verified_rows=len(rows),missing=0,duplicates=0)


def protected_output(source,out):
    code=Path(__file__).resolve().parent
    for parent in (source,code):
        if out==parent or out.is_relative_to(parent) or parent.is_relative_to(out):
            raise ValueError('Report output must not overlap source results or the code directory')


def verify_finished(out,repair_id):
    manifest=read_json(out/'repair_manifest.json')
    if manifest.get('repair_id')!=repair_id or manifest.get('status')!='COMPLETED':
        raise ValueError('Output already exists with different/incomplete content; choose another --output')
    for name,record in manifest.get('outputs',{}).items():
        p=(out/name).resolve()
        if not p.is_relative_to(out) or not p.is_file() or p.stat().st_size!=record['bytes'] or sha(p)!=record['sha256']:
            raise ValueError('Previously repaired report was modified: '+name)
    if not manifest.get('outputs'):raise ValueError('Empty report manifest')
    return manifest


def package(out):
    target=out.with_name(out.name+'.zip')
    files={str(p.relative_to(out)).replace('\\','/'):p for p in out.rglob('*') if p.is_file()}
    if target.exists():
        with zipfile.ZipFile(target) as z:
            if len(z.namelist())!=len(files) or set(z.namelist())!=set(files):raise ValueError('Existing ZIP inventory differs')
            import hashlib
            for name,p in files.items():
                if hashlib.sha256(z.read(name)).hexdigest()!=sha(p):raise ValueError('Existing ZIP was modified')
        return target
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,p in sorted(files.items()):z.write(p,arcname=name)
    return target


def repair(source,output=None,include_dataset_metrics=False):
    source=Path(source).expanduser().resolve(strict=True)
    with input_read_lock(source):
        protocol=load_source(source)
        paths=[source/'protocol.json',source/'state.json']
        for name in ('environment.json','calibration_info.json'):
            if (source/name).is_file():paths.append(source/name)
        supplied={}
        for phase in ('gaussian','counts'):
            if any(j['phase']==phase for j in protocol['jobs']):
                path=source/'review'/f'{phase}_dataset_metrics.csv'
                if not path.is_file():raise ValueError('Missing exported dataset metrics: '+str(path))
                paths.append(path)
        optional=['analytic_sensitivity.csv','heldout_rule_evaluation.csv','diagnostics.csv',
                  'retry_attempts.csv','messages.csv','unavailable_genes.csv','failure_counts.csv']
        for name in optional:
            path=source/'review'/name
            if path.is_file():paths.append(path)
        initial=snapshot(source,paths)
        code_root=Path(__file__).resolve().parent
        code={str(p.relative_to(code_root)).replace('\\','/'):sha(p) for p in sorted(code_root.rglob('*'))
              if p.suffix in ('.py','.R') and '__pycache__' not in p.parts}
        identity_info=dict(source_run_id=protocol['run_id'],inputs=initial,code=code,
                           numpy=np.__version__,scipy=scipy.__version__,include_dataset_metrics=include_dataset_metrics)
        repair_id=identity(identity_info)
        out=(Path(output).expanduser() if output else source.parent/(source.name+'_report_fixed_'+repair_id[:10])).resolve()
        protected_output(source,out)
        if out.exists():
            verify_finished(out,repair_id)
            print('REUSE VERIFIED REPORT',out,flush=True)
            print('TAKE ONLY',package(out),flush=True)
            return out
        checks={}
        genes=read_json(source/'calibration_info.json').get('used_genes')
        for phase in ('gaussian','counts'):
            if not any(j['phase']==phase for j in protocol['jobs']):continue
            print('VERIFY',phase,'exported metrics',flush=True)
            rows=csv_rows(source/'review'/f'{phase}_dataset_metrics.csv')
            checks[phase]=validate_export(rows,protocol['jobs'],phase,protocol,genes)
            supplied[phase]=rows
        sensitivity=csv_rows(source/'review/analytic_sensitivity.csv') if (source/'review/analytic_sensitivity.csv').is_file() else []
        out.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.p1_report_repair_',dir=out.parent) as temp:
            stage=Path(temp)/'result';(stage/'review').mkdir(parents=True)
            issues=build_report(stage,protocol,supplied_rows=supplied,sensitivity_rows=sensitivity,
                                include_dataset_metrics=include_dataset_metrics)
            diagnostics_available=[]
            for name in optional[2:]:
                path=source/'review'/name
                if path.is_file():
                    write_csv(stage/'review'/name,csv_rows(path));diagnostics_available.append(name)
            old_path=source/'review/heldout_rule_evaluation.csv'
            old=csv_rows(old_path) if old_path.is_file() else []
            new=csv_rows(stage/'review/heldout_rule_evaluation.csv') if (stage/'review/heldout_rule_evaluation.csv').is_file() else []
            lookup={(r['design'],int(r['scenario'])):r for r in old}
            changes=[]
            for r in new:
                previous=lookup.get((r['design'],int(r['scenario'])))
                if previous is not None and (previous.get('model')!=r.get('model') or previous.get('status')!=r.get('status')):
                    changes.append(dict(design=r['design'],scenario=r['scenario'],old_model=previous.get('model'),
                                        new_model=r.get('model'),old_status=previous.get('status'),new_status=r.get('status')))
            write_csv(stage/'review/selection_changes.csv',changes)
            checks.update(source_completed=True,source_signature_verified=True,model_fits_rerun=0,
                          historical_selected=sum(bool(r.get('model')) for r in old) if old else None,
                          corrected_selected=sum(bool(r.get('model')) for r in new),
                          corrected_status_counts=dict(Counter(r['status'] for r in new)),
                          changes=changes,whole_arm_failures=issues,
                          diagnostic_tables_available=diagnostics_available,
                          historical_diagnostics_reconstructed=False)
            atomic_json(stage/'checks.json',checks)
            atomic_json(stage/'source_protocol.json',protocol)
            lines=['# Paper 1 report repair','',
                   'No models were refitted. Source results and frozen source were not modified.',
                   f'Source run: {protocol["run_id"]}',f'Source version: {protocol["version"]}',
                   f'Report code version: {VERSION}',
                   f'Historical conditional selections: {checks["historical_selected"]}',
                   f'Corrected conditional selections: {checks["corrected_selected"]}',
                   '','## Statistical changes',
                   '- Global-null FDP uses exact binomial uncertainty for every event count.',
                   '- Non-null fractional FDP retains approximate t intervals; no binomial substitution.',
                   '- Held-out decisions show uncertainty and nominal/screening threshold exceedances.',
                   '- The selector still requires generating context; unknown context means abstention.',
                   '','## Data and diagnostics',
                   '- Dataset/model keys and reported arithmetic were checked against the original protocol.',
                   '- Missing historical diagnostic records cannot be reconstructed from aggregate metrics.',
                   '- Dataset-level exports remain in the source package unless explicitly included here.',
                   '- checks.json and repair_manifest.json identify the inputs and generated outputs.']
            (stage/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
            if snapshot(source,paths)!=initial:raise ValueError('Source changed during report repair; nothing published')
            manifest=dict(status='COMPLETED',kind='report_only_repair',repair_id=repair_id,
                          source_location=str(source),source_results_modified=False,**identity_info)
            manifest['outputs']=snapshot(stage,[p for p in stage.rglob('*') if p.is_file()])
            atomic_json(stage/'repair_manifest.json',manifest)
            if out.exists():raise ValueError('Output appeared during repair; refusing to overwrite it')
            stage.rename(out)
        print('REPORT',out/'review/report.html',flush=True)
        print('CONDITIONAL SELECTIONS',checks['historical_selected'],'->',checks['corrected_selected'],flush=True)
        print('MODEL FITS RERUN 0',flush=True)
        print('TAKE ONLY',package(out),flush=True)
        return out


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-results',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--include-dataset-metrics',action='store_true')
    a=parser.parse_args(argv)
    repair(a.previous_results,a.output,a.include_dataset_metrics)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
