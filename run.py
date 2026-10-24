#!/usr/bin/env python3
"""Reproduce the Paper1 expression-model comparison and validation analyses."""
import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task', choices=['all', 'prepare', 'core', 'matched', 'benchmark', 'supplementary',
                                     'figures', 'report', 'check', 'install-r'], default='all')
    p.add_argument('--primary-only', action='store_true', help='Omit supplementary analyses explicitly')
    p.add_argument('--output', type=Path, default=ROOT / 'results')
    p.add_argument('--raw-cache', type=Path, action='append', default=[],
                   help='Existing download cache; never imports prepared data or results')
    p.add_argument('--resume', action='store_true', help='Explicitly continue this same raw-input experiment')
    p.add_argument('--workers', type=int, default=min(36, os.cpu_count() or 1))
    p.add_argument('--rscript', default=shutil.which('Rscript') or 'Rscript')
    p.add_argument('--r-library', type=Path, default=ROOT / 'R_library')
    p.add_argument('--setup', action='store_true', help='Install missing R dependencies before analysis')
    p.add_argument('--plan', action='store_true', help='Print commands without running analyses')
    p.add_argument('--status', action='store_true', help='Read-only progress of all analysis components')
    p.add_argument('--coverage', action='store_true', help='Report implemented and missing manuscript analyses')
    p.add_argument('--require-complete-paper', action='store_true',
                   help='Stop before computation if any manuscript analysis lacks an integrated workflow')
    return p


def commands(args):
    out = args.output.expanduser().resolve()
    library = args.r_library.expanduser().resolve()
    exe = shutil.which(args.rscript) or args.rscript
    core = [sys.executable, '-m', 'pipeline.core', '--work', str(out / 'core'),
            '--library', str(library), '--rscript', str(exe), '--workers', str(min(6, args.workers))]
    core_started = args.resume and (out / 'core' / 'fresh_run.json').is_file()
    if core_started:
        core += ['--resume']
    for cache in args.raw_cache:
        core += ['--raw-cache', str(cache.expanduser().resolve())]
    jobs = []
    if args.setup or args.task == 'install-r':
        jobs.append(core + ['--task', 'install-r'])
        if '--resume' not in core:
            core += ['--resume']
    if args.task in ('all', 'prepare', 'core', 'check'):
        task = {'all': 'all', 'prepare': 'prepare', 'core': 'models', 'check': 'check'}[args.task]
        cmd = core + ['--task', task]
        jobs.append(cmd)
    if args.task in ('all', 'matched'):
        source = out / 'core'
        jobs.append([sys.executable, '-m', 'pipeline.matched', '--data-work', str(source.resolve()),
                     '--output', str(out / 'matched'), '--workers', str(args.workers),
                     '--rscript', str(exe), '--r-library', str(library), '--offline'] +
                    (['--resume'] if args.resume and (out / 'matched' / 'fresh_run.json').is_file() else []))
    if args.task in ('all', 'benchmark'):
        jobs.append([sys.executable, '-m', 'pipeline.benchmark', '--output', str(out / 'benchmark'),
                     '--data-work', str(out / 'core'),
                     '--workers', str(args.workers), '--r-workers', str(args.workers),
                     '--rscript', str(exe), '--r-library', str(library), '--offline'] +
                    (['--resume'] if args.resume and (out / 'benchmark' / 'fresh_run.json').is_file() else []))
    if args.task == 'supplementary' or args.task == 'all' and not args.primary_only:
        supplementary = [sys.executable, '-m', 'pipeline.supplementary.driver', '--core', str(out / 'core'),
            '--output', str(out / 'supplementary'), '--workers', str(args.workers),
            '--rscript', str(exe), '--r-library', str(library)]
        for cache in args.raw_cache:
            supplementary += ['--raw-cache', str(cache.expanduser().resolve())]
        if args.resume and (out / 'supplementary' / 'fresh_run.json').is_file():
            supplementary += ['--resume']
        jobs.append(supplementary)
    if args.task in ('all', 'figures'):
        jobs.append([sys.executable, '-m', 'pipeline.figures', '--core', str(out / 'core'),
                     '--matched', str(out / 'matched'), '--benchmark', str(out / 'benchmark'),
                     '--output', str(out / 'figures')])
    if args.task == 'report' or args.task == 'all' and not args.primary_only:
        jobs.append([sys.executable, '-m', 'pipeline.reporting', '--run', str(out)])
    return jobs


def main(argv=None):
    args = parser().parse_args(argv)
    if args.coverage:
        import json
        from pipeline.coverage import report
        print(json.dumps(report(), ensure_ascii=False, indent=2))
        return 0
    if args.require_complete_paper:
        from pipeline.coverage import require_complete
        require_complete()
    if args.workers < 1:
        raise ValueError('--workers must be positive')
    out = args.output.expanduser().resolve()
    if args.status:
        from pipeline.matched import read_json, status
        from pipeline.benchmark import status as benchmark_status
        print('Overall:', read_json(out / 'state.json'))
        print('Core:', read_json(out / 'core' / 'status.json'))
        print('Supplementary:', read_json(out / 'supplementary' / 'status.json'))
        print('Reporting:', read_json(out / 'reporting' / 'status.json'))
        status(out / 'matched')
        benchmark_status(out / 'benchmark')
        return 0
    from pipeline.configuration import VERSION, load_config
    load_config(ROOT / 'config.json', args.workers)
    jobs = commands(args)
    if args.plan:
        for command in jobs:
            print(shlex.join(command))
        return 0
    if sys.platform != 'linux' and args.task != 'figures':
        raise RuntimeError('Scientific execution requires Linux. --plan and Python tests work on Windows.')
    if args.task == 'supplementary' or args.task == 'all' and not args.primary_only:
        from pipeline.supplementary.driver import check_python, GROUPS
        check_python(GROUPS)
    from pipeline.matched import lock_output, atomic_json, read_json
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    env['PYTHONPATH'] = str(ROOT) + os.pathsep + env.get('PYTHONPATH', '')
    with lock_output(out):
        from pipeline.lifecycle import claim, code_identity, assert_code
        from pipeline.coverage import report
        identity = code_identity(ROOT)
        run_id = claim(out, 'paper1', args.resume, code=identity)
        env['PAPER1_FRESH_RUN_ID'] = run_id
        atomic_json(out / 'owner.json', {'runner': 'paper1', 'version': VERSION})
        coverage = report()
        atomic_json(out / 'paper_coverage.json', coverage)
        state = {'status': 'RUNNING', 'task': args.task, 'workers': args.workers,
                 'pid': os.getpid(), 'run_id': run_id, 'code': identity, 'started': time.time(),
                 'complete_paper_coverage': coverage['complete_paper']}
        atomic_json(out / 'state.json', state)
        try:
            for command in jobs:
                assert_code(ROOT, identity)
                if not ('--task' in command and command[command.index('--task') + 1] == 'install-r'):
                    from pipeline.environment import snapshot, freeze_environment
                    freeze_environment(out / 'software.json', snapshot(args.rscript, [args.r_library]))
                state['component'] = command[2]
                atomic_json(out / 'state.json', state)
                print('RUN', shlex.join(command), flush=True)
                subprocess.run(command, cwd=ROOT, env=env, check=True)
                assert_code(ROOT, identity)
                if not ('--task' in command and command[command.index('--task') + 1] == 'install-r'):
                    freeze_environment(out / 'software.json', snapshot(args.rscript, [args.r_library]))
            state['status'] = (('PRIMARY_ANALYSES_COMPLETE' if args.primary_only else 'INTEGRATED_ANALYSES_COMPLETE') if args.task == 'all'
                               else 'REQUESTED_TASK_COMPLETE')
            state['finished'] = time.time()
            atomic_json(out / 'state.json', state)
        except BaseException as exc:
            state.update(status='INTERRUPTED' if isinstance(exc, KeyboardInterrupt) else 'FAILED',
                         reason=str(exc), finished=time.time())
            atomic_json(out / 'state.json', state)
            raise
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
