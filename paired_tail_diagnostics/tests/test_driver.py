from pathlib import Path
import tempfile
import unittest
from paired_tail_diagnostics import driver as d


class DriverTests(unittest.TestCase):
    def test_grid(self):
        self.assertEqual([s['id'] for s in d.scenarios()],['s01','s09','s11','s12','s15','s23'])
        self.assertEqual(len(d.KEYS),12)

    def test_seeds_separated(self):
        old={20261107+s['index']*100000+r for s in d.b.scenarios('full') for r in range(1,1001)}
        new={120261108+s['index']*100000+r for s in d.scenarios() for r in range(1,501)}
        aux={x+50000000 for x in new};gauss={x+100000000 for x in new}
        self.assertEqual(len(new),3000)
        self.assertTrue(new.isdisjoint(old))
        self.assertEqual(len(new|aux|gauss),9000)

    def test_job_is_not_completion(self):
        with tempfile.TemporaryDirectory() as t:
            p=dict(output=t,run_id='test');s=d.scenarios()[0]
            d.b.atomic_json(d.b.folder(p,s,1)/'job.json',dict(id=1))
            self.assertFalse(d.valid(p,s,1))

    def test_empty_receipt(self):
        with tempfile.TemporaryDirectory() as t:
            p=dict(output=t,run_id='test');s=d.scenarios()[0]
            d.b.atomic_json(d.b.folder(p,s,1)/'receipt.json',dict(run_id='test',files={}))
            self.assertFalse(d.valid(p,s,1))

    def test_source_files_present(self):
        self.assertEqual(len(set(d.FILES)),len(d.FILES))
        self.assertTrue(all((d.ROOT/n).is_file() for n in d.FILES))

if __name__=='__main__':unittest.main()
