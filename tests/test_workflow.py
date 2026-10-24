"""Public entry points, simulation scope and resumable task validation."""
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import run
from pipeline import benchmark, matched


class WorkflowTests(unittest.TestCase):
    def test_full_plan(self):
        args=run.parser().parse_args(['--output','/tmp/paper1','--workers','60','--raw-cache','/tmp/data'])
        jobs=run.commands(args)
        self.assertEqual([j[2] for j in jobs],['pipeline.core','pipeline.matched','pipeline.benchmark',
                                             'pipeline.supplementary.driver','pipeline.figures','pipeline.reporting'])
        self.assertIn('all',jobs[0])
        self.assertIn('--raw-cache',jobs[0])
        self.assertNotIn('--resume',jobs[0])
        self.assertNotIn('--data-work',jobs[0])
        self.assertIn('--data-work',jobs[2])
        self.assertEqual(jobs[1][jobs[1].index('--workers')+1],'60')
        self.assertEqual(jobs[2][jobs[2].index('--r-workers')+1],'60')

    def test_explicit_workers_not_silently_reduced(self):
        with patch.object(matched,'resources',return_value={'cpus':72,'available_GiB':138,'load_1min':33}):
            self.assertEqual(matched.worker_limit(60),60)
            with self.assertRaisesRegex(RuntimeError,'affinity'):
                matched.worker_limit(73)

    def test_count_registry(self):
        jobs=matched.make_jobs({'contexts':[str(i) for i in range(24)],
                                'mouse_specs':[str(i) for i in range(6)]},
                               SimpleNamespace(reps=200,batch_size=2))
        self.assertEqual(len(jobs),830)
        count=[j for j in jobs if j['phase']=='counts']
        self.assertEqual(sum(j['end']-j['start'] for j in count),1600)
        self.assertEqual(len({j['id'] for j in jobs}),len(jobs))

    def test_partial_task_is_not_complete(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);job={'id':'counts_s00_r00000'}
            (out/'tasks'/job['id']).mkdir(parents=True)
            (out/'tasks'/job['id']/'metrics.csv').write_text('partial')
            self.assertFalse(matched.receipt_ok(out,job,'protocol'))

    def test_damaged_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);job={'id':'x'};(out/'receipts').mkdir()
            (out/'a.csv').write_text('abc')
            record={'job':job,'protocol_hash':'p','files':{'a.csv':{'bytes':3,'sha256':'incorrect'}}}
            (out/'receipts/x.json').write_text(json.dumps(record))
            with self.assertRaisesRegex(RuntimeError,'Damaged'):
                matched.receipt_ok(out,job,'p')

    def test_fixed_design_identities(self):
        benchmark.self_test()

    def test_removed_prepared_import_option(self):
        with self.assertRaises(SystemExit):
            run.parser().parse_args(['--data-work','/tmp/old'])

    def test_setup_keeps_one_experiment(self):
        jobs = run.commands(run.parser().parse_args(['--setup', '--output', '/tmp/fresh']))
        self.assertNotIn('--resume', jobs[0])
        self.assertIn('--resume', jobs[1])


if __name__=='__main__':
    unittest.main()
