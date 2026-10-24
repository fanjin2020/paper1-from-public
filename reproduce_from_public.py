#!/usr/bin/env python3
"""Select existing compatible environments and start a new public-input reproduction."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import site
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def distinct_existing(paths, *, resolve_symlinks=True):
    result = []
    for path in paths:
        p = Path(path).expanduser()
        p = p.resolve() if resolve_symlinks else p.absolute()
        if p.exists() and p not in result:
            result.append(p)
    return result


def python_candidates(explicit=None):
    # A venv's Python symlink must be invoked in place to retain its environment.
    if explicit:
        return distinct_existing([explicit], resolve_symlinks=False)
    values = [sys.executable]
    for base in (ROOT, ROOT.parent / 'paper1-from-publicV1', ROOT.parent / 'paper1_public_evidence'):
        for layout in (base, base / 'paper1-from-public'):
            for name in ('.raw-reproduction-venv', '.venv', '.benchmark-python'):
                values.extend([layout / name / 'bin/python', layout / name / 'Scripts/python.exe'])
    base = ROOT.parent / 'paper1_public_evidence'
    if base.is_dir():
        for pattern in ('paper1*source*/.venv/bin/python', 'paper3*source*/.venv/bin/python'):
            values.extend(sorted(base.glob(pattern)))
    return distinct_existing(values, resolve_symlinks=False)


def optional_package_paths(explicit=()):
    if explicit:
        return distinct_existing(explicit)
    values = list(site.getsitepackages())
    related = ROOT.parent / 'paper1_public_evidence'
    if related.is_dir():
        for pattern in ('paper1*source*/.venv/lib/python3.12/site-packages',
                        'paper3*source*/.venv/lib/python3.12/site-packages'):
            values.extend(sorted(related.glob(pattern)))
    return distinct_existing(values)


def select_python(explicit=None, supplementary=False, package_paths=()):
    code = ('import sys,json; sys.path.insert(0,sys.argv[1]); '
            'sys.path.extend(json.loads(sys.argv[2])); '
            'from pipeline.environment import python_snapshot; '
            'record=python_snapshot(); ')
    if supplementary:
        code += ('from pipeline.supplementary.driver import check_python,GROUPS; '
                 'record["supplementary_packages"]=check_python(GROUPS); ')
    code += 'print("PAPER1_PYTHON\\t"+json.dumps(record))'
    for exe in python_candidates(explicit):
        try:
            result = subprocess.run([str(exe), '-B', '-c', code, str(ROOT), json.dumps(list(map(str, package_paths)))], capture_output=True,
                                    text=True, encoding='utf-8', errors='replace', timeout=180)
        except (OSError, subprocess.TimeoutExpired) as exc:
            print('PYTHON NOT SELECTED', exe, str(exc), flush=True)
            continue
        lines = [line for line in result.stdout.splitlines() if line.startswith('PAPER1_PYTHON\t')]
        if result.returncode == 0 and len(lines) == 1:
            try:
                return exe, json.loads(lines[0].split('\t', 1)[1])
            except ValueError:
                print('PYTHON NOT SELECTED', exe, 'Invalid probe record', flush=True)
                continue
        print('PYTHON NOT SELECTED', exe, result.stderr[-2000:], flush=True)
    raise RuntimeError('No compatible existing Python found. No packages were downloaded. '
                       'Use --python /path/to/an/existing/python matching requirements.txt.')


def library_candidates(explicit=()):
    values = list(explicit)
    # This is a library search only. No results or prepared data enter the analysis.
    related = ROOT.parent / 'paper1_public_evidence'
    values += [related / 'paper1_two_analyses/R_library', ROOT / 'R_library']
    for base in (ROOT.parent, related):
        if base.is_dir():
            for child in sorted(base.iterdir()):
                if child.is_dir() and child.name.startswith('paper1'):
                    values += [child / 'R_library', child / 'runtime/R_library']
    for name in ('R_LIBS', 'R_LIBS_USER', 'PAPER1_R_LIBRARIES'):
        values += [x for x in os.environ.get(name, '').split(os.pathsep) if x]
    return distinct_existing(values)


def select_r(executable=None, explicit=()):
    from pipeline.environment import r_snapshot
    exe = executable or shutil.which('Rscript')
    if not exe:
        raise RuntimeError('Rscript not found. Activate your existing paper1 environment or use --rscript.')
    libraries = library_candidates(explicit)
    record = r_snapshot(exe, libraries)
    return str(exe), libraries, record


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    p.add_argument('--workers', type=int, default=min(60, len(os.sched_getaffinity(0))
                   if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1))
    p.add_argument('--raw-cache', type=Path, action='append', default=[])
    p.add_argument('--python', dest='python_exe')
    p.add_argument('--python-package-path', type=Path, action='append', default=[],
                   help='Existing optional-package site-packages directory; checked, never installed')
    p.add_argument('--rscript')
    p.add_argument('--r-library', action='append', default=[])
    p.add_argument('--check-env', action='store_true', help='Check only; no data downloads or analyses')
    p.add_argument('--primary-only', action='store_true',
                   help='Explicitly run only the implemented central analyses, not the entire manuscript')
    p.add_argument('--resume', action='store_true', help='Continue this same run; requires --output')
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    from pipeline.coverage import require_complete, report
    if not a.check_env and not a.primary_only:
        require_complete()
    if a.workers < 1 or (a.resume and a.output is None):
        raise ValueError('Use positive workers; --resume requires the exact existing --output.')
    from pipeline.configuration import load_config
    config = load_config(ROOT / 'config.json', a.workers)
    print('CONFIG', config['analysis_version'], ROOT / 'config.json', flush=True)
    print('PACKAGE DOWNLOADS DISABLED', flush=True)
    for path in a.python_package_path:
        if not path.is_dir():
            raise ValueError('Existing package directory not found: ' + str(path))
    package_paths = optional_package_paths(a.python_package_path) if not a.primary_only else []
    python, py_info = select_python(a.python_exe, supplementary=not a.primary_only,
                                   package_paths=package_paths)
    rscript, libraries, r_info = select_r(a.rscript, a.r_library)
    print('PYTHON', python, flush=True)
    print('Rscript', rscript, flush=True)
    for p in libraries:
        print('R_LIBRARY', p, flush=True)
    if a.check_env:
        print(json.dumps({'configuration': {'path': str(ROOT / 'config.json'),
                                           'analysis_version': config['analysis_version']},
                          'environment': {'python': py_info, 'R': r_info},
                          'complete_paper_coverage': report()['complete_paper']}, indent=2))
        return 0
    if sys.platform != 'linux':
        raise RuntimeError('Scientific execution requires Linux; environment checks work on Windows.')
    output = (a.output or ROOT.parent / ('paper1_from_public_' + datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    command = [str(python), '-B', str(ROOT / 'run.py'), '--output', str(output),
               '--workers', str(a.workers), '--rscript', rscript]
    if libraries:
        command += ['--r-library', str(libraries[0])]
    for raw in a.raw_cache:
        command += ['--raw-cache', str(raw.expanduser().resolve())]
    if a.resume:
        command += ['--resume']
    if a.primary_only:
        command += ['--primary-only']
    env = os.environ.copy()
    env['R_LIBS'] = os.pathsep.join(map(str, libraries))
    env['PYTHONPATH'] = str(ROOT)
    # The selected interpreter's own site-packages must retain precedence over additions.
    if package_paths:
        code = ('import sys,runpy,json; sys.path.extend(json.loads(sys.argv.pop(1))); '
                'entry=sys.argv.pop(1); sys.argv[0]=entry; runpy.run_path(entry,run_name="__main__")')
        command = [str(python), '-B', '-c', code, json.dumps(list(map(str, package_paths))), *command[2:]]
        env['PAPER1_OPTIONAL_PACKAGE_PATHS'] = json.dumps(list(map(str, package_paths)))
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    print('SCOPE', 'PRIMARY_ANALYSES_ONLY' if a.primary_only else 'INTEGRATED_PUBLIC_INPUT_ANALYSES', flush=True)
    print('OUTPUT', output, flush=True)
    return subprocess.call(command, cwd=ROOT, env=env)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError) as exc:
        print('STOP:', exc, file=sys.stderr)
        raise SystemExit(2)
