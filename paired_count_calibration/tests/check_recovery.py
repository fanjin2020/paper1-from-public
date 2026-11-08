"""Opt-in real-R integration test. Uses a new temporary directory only."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from paired_count_calibration import driver as d


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rscript',default='Rscript')
    ap.add_argument('--r-library',action='append',default=[])
    args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='p1_paired_recovery_') as temporary:
        out=Path(temporary)/'run'
        base=[sys.executable,str(ROOT/'paper1_paired_count.py'),'--output',str(out),
              '--rscript',args.rscript,'--allow-version-drift']
        for lib in args.r_library:base.extend(['--r-library',lib])
        def call(extra,success=True):
            r=subprocess.run(base+extra,capture_output=True,text=True,encoding='utf-8',errors='replace')
            if success and r.returncode:raise AssertionError(r.stdout+'\n'+r.stderr)
            if not success and r.returncode==0:raise AssertionError('Unexpected success')
            return r.stdout+'\n'+r.stderr
        call(['--profile','smoke','--workers','2'])
        p=d.read(out/'protocol.json')
        assert d.progress(p)[0]==4
        snapshots={s['id']:d.read(d.folder(p,s,1)/'result.json') for s in p['scenarios']}
        hashes={s['id']:d.sha(d.folder(p,s,1)/'result.json') for s in p['scenarios']}
        call(['--resume','--workers','3'])
        assert all(d.sha(d.folder(p,s,1)/'result.json')==hashes[s['id']] for s in p['scenarios'])
        s0,s1,s2,_=p['scenarios']
        r=copy.deepcopy(snapshots[s0['id']]);r['metrics'][0]['fpr_050']=.999
        d.atomic_json(d.folder(p,s0,1)/'result.json',r)
        (d.folder(p,s1,1)/'receipt.json').rename(d.folder(p,s1,1)/'receipt.saved')
        for name in ('receipt.json','result.json'):
            (d.folder(p,s2,1)/name).rename(d.folder(p,s2,1)/(name+'.saved'))
        d.atomic_json(d.folder(p,s2,1)/'job.json',dict(interrupted=True))
        status=call(['--status'])
        assert 'VERIFIED DATASETS 1/4' in status,status
        call(['--resume','--workers','2'])
        assert d.progress(p)[0]==4
        for s in p['scenarios']:
            before=snapshots[s['id']];after=d.read(d.folder(p,s,1)/'result.json')
            before.pop('attempt',None);after.pop('attempt',None)
            assert before==after,'Seeded scientific output changed on resume'
        message=call(['--resume','--genes','251'],success=False)
        assert 'Resume setting mismatch: genes' in message,message
        assert (out/'review.zip').is_file()
    print('PASS: real R pipeline, complete resume, corruption, missing receipt, job-only recovery, seeded equality, settings guard')


if __name__=='__main__':main()
