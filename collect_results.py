#!/usr/bin/env python3
"""Collect small completed-result exports without counts, raw inputs or ZIP archives."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(argv);root=a.run.expanduser().resolve();out=a.output.expanduser().resolve()
    if root==out or out.is_relative_to(root) or root.is_relative_to(out):
        raise ValueError('Collection output must be separate from the experiment directory')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Collection output must be empty')
    state=json.loads((root/'state.json').read_text(encoding='utf-8'))
    report=json.loads((root/'reporting/status.json').read_text(encoding='utf-8'))
    if state.get('status')!='INTEGRATED_ANALYSES_COMPLETE' or report.get('status')!='RESULT_EXPORTS_COMPLETE':
        raise ValueError('The integrated analysis and result exports are not complete')
    if state.get('run_id')!=report.get('run_id'):
        raise ValueError('Result identities differ')
    with (root/'reporting/source_index.csv').open(encoding='utf-8',newline='') as stream:
        for record in csv.DictReader(stream):
            source=(root/'reporting'/record['file']).resolve()
            if not source.is_relative_to(root/'reporting') or not source.is_file():
                raise ValueError('Invalid export path: '+record['file'])
            digest=hashlib.sha256()
            with source.open('rb') as data:
                for block in iter(lambda:data.read(8*1024**2),b''):digest.update(block)
            if source.stat().st_size!=int(record['bytes']) or digest.hexdigest()!=record['sha256']:
                raise ValueError('Export changed after completion: '+record['file'])
    out.mkdir(parents=True,exist_ok=True)
    for name in ['reporting','figures']:
        shutil.copytree(root/name,out/name)
    for name in ['state.json','software.json','fresh_run.json','paper_coverage.json']:
        shutil.copyfile(root/name,out/name)
    for folder,name in [('matched/results','analysis_summary.json'),('supplementary','python.json'),
                        ('supplementary','status.json'),('benchmark','run_state.json')]:
        source=root/folder/name
        if not source.is_file():raise ValueError('Missing completion record: '+str(source))
        target=out/folder/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
    total=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    print('RESULTS',out)
    print('SIZE_MiB',round(total/1024**2,1))
    print('No raw counts, prepared matrices, cell-level annotations, or simulation batches copied.')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError) as exc:
        print('STOP:',exc,file=sys.stderr);raise SystemExit(2)
