import csv
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from publication_calibration import report


class CalibrationPublicationTests(unittest.TestCase):
    def test_complete_inputs(self):
        self.assertEqual(len(report.verify()), 6)
        self.assertEqual(len(report.summarize()), 16)

    def test_gaussian_headline(self):
        rows = report.read('tail_summary')
        r = next(r for r in rows if (r['scenario'], r['arm'], r['test']) == ('s09', 'gaussian_unit', 'HC3'))
        self.assertAlmostEqual(float(r['fdp_planned_050']), .108)
        self.assertAlmostEqual(float(r['bonferroni_fwer']), .106)
        self.assertGreater(float(r['coverage95']), .95)
        self.assertGreater(float(r['fdp_planned_050_ci_low']), .05)

    def test_checksums_reject_changed_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'source'
            shutil.copytree(report.ROOT, root)
            path = root / 'data/count_summary.csv.gz'
            path.write_bytes(path.read_bytes() + b'x')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                report.verify(root)

    def test_output_safety_and_separate_families(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'report'
            self.assertEqual(report.main(['--output', str(out)]), 0)
            rows = json.loads((out / 'summary.json').read_text())
            r = next(r for r in rows if (r['study'], r['arm'], r['test']) == ('tail', 'gaussian_unit', 'HC3'))
            self.assertGreater(r['fdp_retained_050_max'], r['fdp_planned_050_max'])
            with self.assertRaises(SystemExit):
                report.main(['--output', str(out)])


if __name__ == '__main__':
    unittest.main()
