"""The published held-out decisions must follow the training-only rule."""
import csv
import gzip
from pathlib import Path
import unittest
from design_extension.report import choose, fit_rule

ROOT=Path(__file__).resolve().parents[2]/'publication_figures/data'


def rows(name):
    def typed(v):
        if v=='':return None
        if v in ('True','False'):return v=='True'
        try:return float(v)
        except ValueError:return v
    with gzip.open(ROOT/('extension_'+name+'.csv.gz'),'rt',encoding='utf-8',newline='') as f:
        return [{k:typed(v) for k,v in r.items()} for r in csv.DictReader(f)]


class PublicationEvidenceTests(unittest.TestCase):
    def test_decisions_reproduced_without_holdout_outcomes(self):
        designs=rows('design_diagnostics')
        counts=rows('counts_summary')
        recorded=rows('heldout_rule_evaluation')
        rule=fit_rule(counts,designs)
        lookup={d['design']:d for d in designs}
        self.assertEqual(len(recorded),32)
        for r in recorded:
            actual=choose(rule,lookup[r['design']],r['scenario'])
            self.assertEqual((actual['model'],actual['status']),(r['model'],r['status']))
        self.assertTrue(all(r['split']=='training' for r in rule['training_metrics']))

    def test_unknown_context_remains_an_abstention(self):
        designs=rows('design_diagnostics')
        rule=fit_rule(rows('counts_summary'),designs)
        for d in designs:
            self.assertEqual(choose(rule,d,None)['status'],'generating_context_required')


if __name__=='__main__':unittest.main()
