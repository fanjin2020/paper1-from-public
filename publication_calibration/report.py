"""Summarize completed experiments without refitting or selecting a best test."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parent
ARMS = ('estimated_shared', 'estimated_design', 'independent_shared',
        'independent_design', 'fixed_unit', 'gaussian_unit')
TESTS = ('model_based', 'HC3')


def read(name, root=ROOT):
    with gzip.open(root / 'data' / (name + '.csv.gz'), 'rt', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def verify(root=ROOT):
    manifest = json.loads((root / 'manifest.json').read_text(encoding='utf-8'))
    for item in manifest:
        data = (root / item['file']).read_bytes()
        if hashlib.sha256(data).hexdigest() != item['sha256']:
            raise ValueError('Calibration source checksum mismatch: ' + item['name'])
        rows = list(csv.reader(io.StringIO(gzip.decompress(data).decode('utf-8'))))
        if len(rows) != item['rows'] + 1 or any(len(r) != item['columns'] for r in rows):
            raise ValueError('Calibration source dimensions differ: ' + item['name'])
    for name, dimension, expected, reps in [
        ('count_summary', 'weight_mode', 96, 1000), ('tail_summary', 'arm', 72, 500)
    ]:
        rows = read(name, root)
        keys = {(r['scenario'], r[dimension], r['test']) for r in rows}
        if len(rows) != expected or len(keys) != expected:
            raise ValueError('Duplicate or incomplete grid: ' + name)
        for r in rows:
            if int(r['datasets']) != reps:
                raise ValueError('Unexpected replicate count: ' + name)
            for metric in ('coverage95', 'fpr_050', 'fdp_planned_050', 'fdp_retained_050'):
                value = float(r[metric])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError('Invalid rate: ' + metric)
            if not math.isclose(float(r['coverage95']) + float(r['fpr_050']), 1, abs_tol=1e-10):
                raise ValueError('Coverage and matching two-sided rejection disagree')
    tail = read('tail_summary', root)
    cases = {r['scenario'] for r in tail}
    if cases != {'s01', 's09', 's11', 's12', 's15', 's23'}:
        raise ValueError('Unexpected targeted scenarios')
    if {(r['scenario'], r['arm'], r['test']) for r in tail} != {
        (s, a, t) for s in cases for a in ARMS for t in TESTS
    }:
        raise ValueError('Tail experiment arm grid differs')
    return manifest


def summarize(root=ROOT):
    verify(root)
    result = []
    for study, name, dimension in [('count', 'count_summary', 'weight_mode'),
                                  ('tail', 'tail_summary', 'arm')]:
        rows = read(name, root)
        for arm, test in sorted({(r[dimension], r['test']) for r in rows}):
            group = [r for r in rows if r[dimension] == arm and r['test'] == test]
            item = dict(study=study, arm=arm, test=test, scenarios=len(group))
            for metric in ('coverage95', 'fpr_050', 'fdp_planned_050', 'fdp_retained_050',
                           'bonferroni_fwer'):
                if metric in group[0]:
                    values = [float(r[metric]) for r in group]
                    item[metric + '_min'] = min(values)
                    item[metric + '_max'] = max(values)
            item['pointwise_mc_lower_above_nominal'] = sum(
                float(r['fdp_planned_050_ci_low']) > .05 for r in group)
            result.append(item)
    return result


def render_markdown(items):
    lines = ['# Paired coefficient test calibration', '',
             'All paired differences are null. FDR is the mean dataset-level false discovery',
             'proportion and equals the probability of any rejection in these experiments.',
             'Ranges below describe scenarios, not confidence intervals. Planned families',
             'contain 3,000 genes; retained families use the expression-filtered set.', '',
             '| Study | Arm | Test | Coverage (%) | Planned BH FDR (%) | Retained BH FDR (%) |',
             '|---|---|---|---:|---:|---:|']
    for r in items:
        def span(k):
            return f"{r[k + '_min'] * 100:.2f}-{r[k + '_max'] * 100:.2f}"
        lines.append(f"| {r['study']} | {r['arm']} | {r['test']} | "
                     f"{span('coverage95')} | {span('fdp_planned_050')} | {span('fdp_retained_050')} |")
    lines += ['', '## Interpretation', '',
              'The current paired HC3 statistic uses the full-model residual t reference.',
              'Coverage near 95% does not establish calibration of its extreme p-value tail.',
              'Inflation persists with unit weights and independent Gaussian responses, and',
              'also occurs with Bonferroni correction. Independent count-derived weights',
              'reduce some discrepancies without universally restoring nominal error rates.',
              'Model-based results are more stable in these settings; weighted count results',
              'still show departures. Neither estimator is automatically selected by this report.', '',
              'Pointwise Monte Carlo intervals quantify simulation uncertainty, not a',
              'simultaneous certification across scenarios. No equivalence or universal',
              'calibration claim is inferred from a nonsignificant excess.', '',
              'The six targeted scenarios were selected after examining the count grid.',
              'Fresh seeds supply new draws, not independent clinical validation.',
              'The Gaussian variance moment-equivalent degrees of freedom describe the',
              'random denominator; they are not substituted into a corrected test.', '',
              'Sources: Tables S22-S23; figures 10, S14 and S15.']
    return '\n'.join(lines) + '\n'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args(argv)
    manifest = verify()
    if args.verify_only:
        print('Verified 24,000 count-grid datasets and 3,000 targeted datasets through compact summaries.')
        return 0
    result = summarize()
    if args.output is None:
        out = Path(tempfile.mkdtemp(prefix='p1_calibration_report_' + datetime.now().strftime('%Y%m%d_%H%M%S_'),
                                    dir=ROOT.parent.parent))
    else:
        out = args.output.resolve()
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            parser.error('Output must be a new or empty directory; source results are read-only')
        out.mkdir(parents=True, exist_ok=True)
    (out / 'calibration_report.md').write_text(render_markdown(result), encoding='utf-8')
    (out / 'summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    for name in ('count_summary', 'tail_summary', 'tail_rates', 'paired_comparisons',
                 'expression_strata', 'gaussian_geometry'):
        (out / (name + '.csv')).write_bytes(gzip.decompress((ROOT / 'data' / (name + '.csv.gz')).read_bytes()))
    (out / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print('REPORT', out)
    print('No downloads, package installs, simulations, refits or adjusted p-value changes.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
