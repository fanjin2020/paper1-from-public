"""Verify compact publication evidence; no downloads or model fitting."""
import csv
import gzip
import hashlib
from pathlib import Path
from publication_tables.export_tables import verify
from publication_calibration.report import verify as verify_calibration

ROOT=Path(__file__).resolve().parent


def rows(name):
    with gzip.open(ROOT/'publication_figures/data'/name,'rt',encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))


def main():
    tables=verify()
    verify_calibration()
    with (ROOT/'publication_figures/source_index.csv').open(encoding='utf-8',newline='') as f:
        inputs=list(csv.DictReader(f))
    for row in inputs:
        path=ROOT/'publication_figures'/row['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Figure input checksum mismatch: '+row['file'])
    held=rows('extension_heldout_rule_evaluation.csv.gz')
    statuses={s:sum(r['status']==s for r in held) for s in
              ('conditionally_eligible','no_supported_selection','within_effect_not_estimable')}
    expected={'conditionally_eligible':16,'no_supported_selection':8,'within_effect_not_estimable':8}
    if statuses!=expected or len({(r['design'],r['scenario']) for r in held})!=32:
        raise ValueError('Held-out combination counts differ')
    selected=[r for r in held if r['status']=='conditionally_eligible']
    if sum(float(r['heldout_fdr'])>.05 for r in selected)!=8:
        raise ValueError('Held-out nominal exceedance count differs')
    if sum(float(r['heldout_fdr_mc_upper'])>.075 for r in selected)!=2:
        raise ValueError('Held-out Monte Carlo screen exceedance count differs')
    for filename,n in [('extension_counts_summary.csv.gz',480),('extension_gaussian_summary.csv.gz',384)]:
        data=rows(filename)
        if len(data)!=n or len({(r['design'],r['scenario'],r['model']) for r in data})!=n:
            raise ValueError('Design benchmark grid is incomplete: '+filename)
    print('Verified',len(tables),'table sections and',len(inputs),'figure inputs')
    print('Held-out: 16 eligible, 8 unsupported, 8 structurally non-estimable')
    print('Figures: 10 main and 15 supplementary; no model fitting performed')


if __name__=='__main__':main()
