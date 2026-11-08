#!/usr/bin/env python3
"""Select an existing numerical environment; never install packages or download inputs."""
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'


def main():
    try:
        import numpy, scipy
    except ImportError:
        candidates=[ROOT/'.raw-reproduction-venv/bin/python',ROOT/'.venv/bin/python']
        base=ROOT.parent/'paper1_public_evidence'
        if base.is_dir():candidates+=sorted(base.glob('paper*source*/.venv/bin/python'))
        for exe in candidates:
            if exe.is_file():
                p=subprocess.run([str(exe),'-c','import numpy,scipy'],capture_output=True)
                if p.returncode==0:
                    print('REUSE PYTHON',exe,flush=True)
                    os.execv(str(exe),[str(exe),str(Path(__file__).absolute()),*sys.argv[1:]])
        raise RuntimeError('numpy/scipy not available. Activate an existing compatible environment; no installs were attempted')
    from design_extension.driver import main as run
    print('NO DOWNLOADS; NO PACKAGE INSTALLS; ORIGINAL RESULTS READ ONLY',flush=True)
    return run()


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,RuntimeError) as e:print('STOP:',e,file=sys.stderr);raise SystemExit(2)
