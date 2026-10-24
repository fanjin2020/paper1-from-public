"""Regression checks for task health and deterministic cancellation directions."""
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pipeline import benchmark as B, matched as M


class TaskHealthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        B.load_science()

    def test_near_zero_tie_is_stable(self):
        for c in ([-.3, 1e-17, 2e-17], [-.3, 2e-17, 1e-17], [-.3, -1e-15, 0.]):
            index, v, tol = B.cancellation_direction(c)
            self.assertEqual(index, 1)
            self.assertLess(abs(np.asarray(c) @ v), 1e-12)
            self.assertGreater(tol, 1e-15)

    def test_non_tie_keeps_smallest_coordinate(self):
        self.assertEqual(B.cancellation_direction([.3, .02, .01])[0], 2)
        with self.assertRaises(ValueError):
            B.cancellation_direction([0., 0., 0.])

    def test_grid_cancellation_truth_and_recorded_direction(self):
        args = B.parser().parse_args(["--module", "paired", "--paired-reps", "12"])
        jobs, designs, _ = B.make_jobs(args)
        self.assertEqual(len(designs), 36)
        for spec in designs:
            g = B.generic_design(spec)
            if g["unit"] < 1e-14:
                continue
            profiles = B.means(g)
            for _, mu, truth in profiles:
                self.assertAlmostEqual(g["d"] @ mu, truth, places=12)
            cancellation = profiles[1][1] - profiles[0][1]
            self.assertAlmostEqual(np.sqrt(np.mean(cancellation ** 2)), .5, places=12)
            self.assertIn("cancellation_coordinate", spec)

    def test_joint_f_batch_invariance(self):
        spec = dict(design="d008", design_index=8, n=40, partial_correlation=.7,
                    leverage="ordinary", weights="uniform", error_law="heteroskedastic_gaussian")
        full = B.paired_batch(dict(spec=spec, start=0, end=40))
        parts = [B.paired_batch(dict(spec=spec, start=a, end=b)) for a, b in ((0, 13), (13, 40))]
        for column in ("joint_F_rejected", "covered", "rejected"):
            np.testing.assert_array_equal(full[column], parts[0][column] + parts[1][column])

    def mouse_fixture(self, folder, failed=None):
        models = ["baseline", "cage_fixed", "cage_correlation", "dream_shared", "dream_native",
                  "mundlak_shared", "mundlak_native"]
        with gzip.open(folder / "background_fits.csv.gz", "wt", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["model", "gene_symbol", "status", "p"])
            writer.writeheader()
            for model in models:
                writer.writerow(dict(model=model, gene_symbol="g", status="fit_failed" if model == failed else "available",
                                     p="" if model == failed else ".5"))

    def test_whole_mouse_arm_rejected_even_with_valid_receipt(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            job = {"id": "mouse_cortex_linear", "phase": "mouse", "spec": "cortex_linear"}
            folder = root / "tasks" / job["id"]
            folder.mkdir(parents=True)
            self.mouse_fixture(folder, "dream_shared")
            file = folder / "background_fits.csv.gz"
            M.atomic_json(root / "receipts" / (job["id"] + ".json"), {
                "job": job, "protocol_hash": "p", "files": {
                    file.relative_to(root).as_posix(): {"sha256": M.sha(file), "bytes": file.stat().st_size}}})
            self.assertTrue(M.receipt_ok(root, job, "p", validate_health=False))
            with self.assertRaisesRegex(RuntimeError, "Whole empirical mouse"):
                M.receipt_ok(root, job, "p")
            self.mouse_fixture(folder)
            M.task_health(folder, job)

    def test_failed_count_model_not_receipted(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            (folder / "metrics.csv").write_text("scenario,replicate,model,usable_run\n" +
                "\n".join(f"s,0,m{i},{'FALSE' if i == 0 else 'TRUE'}" for i in range(5)))
            with self.assertRaisesRegex(RuntimeError, "Whole count model"):
                M.task_health(folder, {"id": "counts_s00_r00000", "phase": "counts", "start": 0, "end": 1})

    def test_summary_flags_not_completion(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            M.atomic_json(root / "results/analysis_summary.json", {"status": "WHOLE_MODEL_FAILURES_REQUIRE_REVIEW"})
            with self.assertRaisesRegex(RuntimeError, "requires review"):
                M.summary_ok(root)
            M.atomic_json(root / "results/analysis_summary.json", {"status": "ANALYSES_COMPLETE"})
            M.summary_ok(root)

    def test_receipt_escape_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "output"
            outside = Path(d) / "outside.csv"
            outside.write_text("x")
            job = {"id": "bad"}
            M.atomic_json(root / "receipts/bad.json", {"job": job, "protocol_hash": "p", "files": {
                "../outside.csv": {"sha256": M.sha(outside), "bytes": 1}}})
            with self.assertRaisesRegex(RuntimeError, "escapes"):
                M.receipt_ok(root, job, "p")

    def test_matched_driver_requests_serial_inner_backend(self):
        source = Path(M.__file__).read_text(encoding="utf-8")
        self.assertIn("inner = 1", source)
        self.assertNotIn("inner = max(1, min(6", source)


if __name__ == "__main__":
    unittest.main()
