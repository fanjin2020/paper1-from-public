#!/usr/bin/env python3
"""Default: repair a completed report. New model fitting requires --run-new."""
import argparse
import os
from pathlib import Path
import sys


DEFAULT_RESULTS=Path('/home/research/git/repository/p1_design_1101_202337')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--run-new',action='store_true',help='Explicitly start a separate new simulation run')
    mode.add_argument('--resume',action='store_true',help='Resume only runs created by this fixed code')
    mode.add_argument('--status',action='store_true',help='Inspect only a run created by this fixed code')
    parser.add_argument('--previous-results',type=Path,default=DEFAULT_RESULTS)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--workers',type=int,default=36)
    parser.add_argument('--include-dataset-metrics',action='store_true')
    args,extra=parser.parse_known_args(argv)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key]='1'
    try:
        import numpy, scipy
    except ImportError as e:
        raise ValueError('Activate the existing paper1 environment with numpy/scipy; no installs were attempted') from e
    if not (args.run_new or args.resume or args.status):
        if extra:parser.error('Unsupported report-repair arguments: '+' '.join(extra))
        from design_extension.repair_report import repair
        print('REPORT REPAIR ONLY; NO MODEL FITTING; NO DOWNLOADS; NO PACKAGE INSTALLS',flush=True)
        repair(args.previous_results,args.output,args.include_dataset_metrics)
        return 0
    if args.output is None:parser.error('--output is required for a new run, resume, or status')
    if args.include_dataset_metrics:parser.error('--include-dataset-metrics applies only to report repair')
    from design_extension import VERSION
    from design_extension.driver import main as run, read_json
    out=args.output.expanduser().resolve()
    if args.status or args.resume:
        p=read_json(out/'protocol.json')
        if not p or p.get('version')!=VERSION:
            raise ValueError('Use the original source for legacy run status/resume. Default mode repairs legacy reports without fitting')
        if extra:parser.error('Resume restores the recorded settings; only --workers may be changed')
        if args.status:return run(['--status','--output',str(out)])
        command=['--resume','--output',str(out),'--workers',str(args.workers)]
        for key in ('profile','phase','count_reps','count_genes','gaussian_reps','seed','fdr_screen'):
            command+=['--'+key.replace('_','-'),str(p[key])]
        if p.get('data_work'):command+=['--data-work',p['data_work']]
        if p.get('rscript'):command+=['--rscript',p['rscript']]
        libs=set()
        for info in p.get('runtime',{}).get('R',{}).get('packages',{}).values():
            if isinstance(info,dict) and info.get('path'):libs.add(str(Path(info['path']).parent))
        for path in sorted(libs):command+=['--r-library',path]
        return run(command)
    if out.exists() and any(out.iterdir()):raise ValueError('--run-new requires a new or empty output directory')
    print('NEW EXPERIMENT EXPLICITLY REQUESTED; ONE COUNT DATASET PER CHECKPOINT',flush=True)
    print('NO DOWNLOADS; NO PACKAGE INSTALLS; EXISTING RESULTS ARE NOT MODIFIED',flush=True)
    return run(['--output',str(out),'--workers',str(args.workers),*extra])


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,RuntimeError,OSError,KeyError) as e:
        print('STOP:',e,file=sys.stderr);raise SystemExit(2)
