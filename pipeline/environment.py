"""Check and record the installed scientific environment without installing packages."""
import importlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

from .io import dump, require

ROOT = Path(__file__).resolve().parents[1]


def python_requirements(root=ROOT):
    requirements = {}
    for line in (Path(root) / 'requirements.txt').read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        match = re.fullmatch(r'([A-Za-z0-9_]+)==([^\s]+)', line)
        require(match is not None, 'Expected an exact Python package pin: ' + line)
        requirements[match[1]] = match[2]
    return requirements


def python_snapshot(root=ROOT):
    require(sys.version_info[:2] == (3, 12), 'Python 3.12 is required for this scientific environment')
    packages, errors = {}, []
    for name, expected in python_requirements(root).items():
        try:
            module = importlib.import_module(name)
            version = str(module.__version__)
            packages[name] = version
            if version != expected:
                errors.append(f'{name}: expected {expected}, found {version}')
        except Exception as exc:
            errors.append(f'{name}: {type(exc).__name__}: {exc}')
    require(not errors, 'Python environment mismatch; no packages were installed.\n' + '\n'.join(errors))
    return {'python': platform.python_version(), 'executable': sys.executable,
            'packages': packages, 'platform': platform.platform()}


def scientific_identity(record):
    """Exclude installation locations, but retain every recorded scientific version."""
    r = record['R']
    return {'python': record['python']['python'], 'python_packages': record['python']['packages'],
            'R': r['R'], 'R_packages': r['packages'],
            'R_loaded_namespaces': r.get('loaded_namespaces', {})}


def freeze_environment(path, record):
    path = Path(path)
    if path.exists():
        previous = json.loads(path.read_text(encoding='utf-8'))
        require(scientific_identity(previous) == scientific_identity(record),
                'Software versions changed within this experiment; start a NEW output directory.')
    else:
        dump(path, record)
    return scientific_identity(record)


def r_snapshot(executable, libraries, root=ROOT):
    root = Path(root)
    result = subprocess.run([str(executable), '--vanilla', str(root / 'analysis/environment.R'),
                             str(root / 'inputs/environment.json'), *map(str, libraries)],
                            capture_output=True, text=True, encoding='utf-8', errors='replace',
                            timeout=300)
    require(result.returncode == 0, 'R environment mismatch; no packages were installed.\n' +
            result.stdout[-6000:] + result.stderr[-6000:])
    records = [line.removeprefix('PAPER1_ENVIRONMENT\t') for line in result.stdout.splitlines()
               if line.startswith('PAPER1_ENVIRONMENT\t')]
    require(len(records) == 1, 'R environment probe did not return exactly one record')
    return json.loads(records[0])


def snapshot(executable, libraries, root=ROOT):
    return {'python': python_snapshot(root), 'R': r_snapshot(executable, libraries, root)}
