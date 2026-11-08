"""Real-R smoke, report and damaged-checkpoint recovery in a temporary directory."""
import argparse
import copy
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from paired_tail_diagnostics import driver as d


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rscript',default='Rscript');ap.add_argument('--r-library',action='append',default=[])
    args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='p1_tail_recovery_') as tmp:
        out=Path(tmp)/'run'
        base=[sys.executable,str(ROOT/'paper1_paired_tail.py'),'--output',str(out),'--rscript',args.rscript,'--allow-version-drift']
        for lib in args.r_library:base+=['--r-library',lib]
        def run(extra,ok=True):
            r=subprocess.run(base+extra,capture_output=True,text=True,encoding='utf-8',errors='replace')
            if ok and r.returncode:raise AssertionError(r.stdout+r.stderr)
            if not ok and not r.returncode:raise AssertionError('Unexpected success')
            return r.stdout+r.stderr
        run(['--profile','pilot','--reps','3','--genes','500','--workers','2'])
        p=d.b.read(out/'protocol.json');assert d.progress(p)[0]==18
        s0,s1,s2=p['scenarios'][:3]
        before={s['id']:d.b.read(d.b.folder(p,s,1)/'result.json') for s in p['scenarios']}
        hashes={s['id']:d.b.sha(d.b.folder(p,s,1)/'result.json') for s in p['scenarios']}
        run(['--resume','--workers','3'])
        assert all(d.b.sha(d.b.folder(p,s,1)/'result.json')==hashes[s['id']] for s in p['scenarios'])
        path=d.b.folder(p,s0,1)/'events.csv.gz'
        with path.open('ab') as f:f.write(b'corrupt')
        (d.b.folder(p,s1,1)/'receipt.json').rename(d.b.folder(p,s1,1)/'receipt.saved')
        for name in ('result.json','events.csv.gz','receipt.json'):
            (d.b.folder(p,s2,1)/name).rename(d.b.folder(p,s2,1)/(name+'.saved'))
        d.b.atomic_json(d.b.folder(p,s2,1)/'job.json',{})
        assert 'VERIFIED 15/18' in run(['--status'])
        run(['--resume','--workers','2']);assert d.progress(p)[0]==18
        for s in p['scenarios']:
            a=before[s['id']];b=d.b.read(d.b.folder(p,s,1)/'result.json')
            a.pop('attempt');b.pop('attempt');assert a==b
        with (out/'summary.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        assert len(rows)==72 and all('bonferroni_fwer' in r for r in rows)
        with (out/'tail_summary.csv').open(newline='') as f:tails=list(csv.DictReader(f))
        assert len(tails)==576
        with (out/'paired_comparisons.csv').open(newline='') as f:pairs=list(csv.DictReader(f))
        assert all((r['arm']=='gaussian_unit')==(r['reference_arm']=='gaussian_unit') for r in pairs)
        assert 'Resume setting mismatch' in run(['--resume','--genes','501'],ok=False)
        assert (out/'review_tail.zip').is_file()
    print('PASS: 18 datasets, 216 method results, reports, unchanged resume, event corruption, missing receipts, job-only recovery, seeded equality')

if __name__=='__main__':main()
