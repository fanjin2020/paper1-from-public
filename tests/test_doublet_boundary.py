import contextlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
from scipy import sparse

from pipeline.supplementary import interaction_math as U, interactions as I


class DoubletBoundaryTests(unittest.TestCase):
    def test_below_minimum_is_unavailable_not_a_cell_decision(self):
        for n in (0, 77, 99):
            for fraction in ('P', 'V'):
                with self.subTest(n=n, fraction=fraction), tempfile.TemporaryDirectory() as d:
                    x = sparse.csr_matrix(np.ones((129, 5), dtype=int))
                    keep = np.arange(129) < n
                    sample = SimpleNamespace(gsm='GSM6911591', donor='donor', fraction=fraction)
                    with patch.dict('sys.modules', {'scrublet': None}):
                        retained, row = U.doublet_risk(x, np.arange(129).astype(str), keep, Path(d)/'lib', sample)
                    self.assertIsNone(retained)
                    self.assertEqual(row['status'], 'not_testable_insufficient_qc_cells')
                    self.assertEqual(row['qc_cells'], n)
                    for field in ('predicted_doublets', 'retained_cells', 'threshold'):
                        self.assertIsNone(row[field])
                    self.assertFalse(row['threshold_review_required'])
                    self.assertEqual(json.loads((Path(d)/'lib/specification.json').read_text()), row)
                    self.assertFalse((Path(d)/'lib/observed_scores.csv.gz').exists())

    def test_invalid_mask_is_still_a_fatal_error(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(U.Blocked, 'mask must be Boolean'):
                U.doublet_risk(sparse.eye(4), ['a','b','c','d'], np.ones(4), Path(d)/'lib', None)

    def test_at_minimum_runs_same_algorithm_and_technical_errors_propagate(self):
        class Scrublet:
            def __init__(self, x, **kwargs):
                assert x.shape[0] == 100
                assert kwargs == {'expected_doublet_rate': .06, 'random_state': 290}
            def scrub_doublets(self, **kwargs):
                assert kwargs == {'use_approx_neighbors': False, 'n_prin_comps': 30, 'verbose': True}
                raise RuntimeError('technical test failure')
        modules = {'scrublet': SimpleNamespace(Scrublet=Scrublet),
                   'threadpoolctl': SimpleNamespace(threadpool_limits=lambda **kwargs: contextlib.nullcontext())}
        with tempfile.TemporaryDirectory() as d, patch.dict('sys.modules', modules):
            with self.assertRaisesRegex(RuntimeError, 'technical test failure'):
                U.doublet_risk(sparse.csr_matrix(np.ones((100, 5))), np.arange(100),
                               np.ones(100, bool), Path(d)/'lib', None)

    def test_missing_filtered_library_is_not_adopted_as_singlet(self):
        # Both matched variants select from the filtered units, not original count columns.
        primary = pd.DataFrame({'donor':['unavailable'], 'unit_id':['u'], 'case':[0]})
        units = pd.DataFrame(columns=['fraction','cell_type','n_cells','donor','unit_id'])
        with self.assertRaises(U.Blocked):
            I.paired_setup(primary, units, 'End')

    def test_driver_records_unavailable_library_and_continues_without_filling_counts(self):
        from pipeline.supplementary import annotation as A, driver as D
        for include_good in (False, True):
            with self.subTest(include_good=include_good), tempfile.TemporaryDirectory() as folder:
                out = Path(folder)
                records = [dict(gsm='small', fraction='V', donor='d0', count_url='small')]
                if include_good:
                    records.append(dict(gsm='good', fraction='P', donor='d1', count_url='good'))
                samples = pd.DataFrame(records)
                units = samples.assign(cell_type='Ast', n_cells=100, unit_id=samples.gsm,
                                       total_umi=300, age=80, sex='F', pmi=3, ph_cerebellum=7, case=0)
                genes = np.array(['g1', 'g2', 'g3'])
                barcodes = np.array([f'cell{i}' for i in range(129)])
                x = sparse.csr_matrix(np.ones((129, 3), dtype=int))
                original = pd.DataFrame(100, index=genes, columns=samples.gsm)
                original.index.name = 'gene_symbol'
                keeps = [np.arange(129) < (77 if s.gsm == 'small' else 100)
                         for s in samples.itertuples()]
                for sample, keep in zip(samples.itertuples(), keeps):
                    pd.DataFrame(dict(barcode=barcodes, qc_pass=keep, label='Ast')).to_csv(
                        out / f'285_{sample.gsm}_cell_annotations.csv.gz', index=False)
                real_risk = U.doublet_risk

                def risk(matrix, bars, keep, directory, sample, annotation):
                    if sample.gsm == 'small':
                        return real_risk(matrix, bars, keep, directory, sample, annotation)
                    return keep, dict(gsm='good', status='completed_risk_screen_not_doublet_truth',
                                      threshold_review_required=False)

                with contextlib.ExitStack() as stack:
                    stack.enter_context(patch.object(A, 'soft', return_value=None))
                    stack.enter_context(patch.object(A, 'target_register', return_value=samples))
                    stack.enter_context(patch.object(A, 'obtain_identity', return_value=None))
                    stack.enter_context(patch.object(A, 'annotate_and_aggregate',
                        return_value=(original, units, 'Ast', None, None)))
                    stack.enter_context(patch.object(A, 'read_counts', return_value=(x, genes, barcodes, None)))
                    stack.enter_context(patch.object(A, 'qc_counts', side_effect=[(k, None, None, None) for k in keeps]))
                    stack.enter_context(patch.object(D, 'panel', return_value=pd.DataFrame(
                        {'gene_symbol': ['g1'], 'discovery_direction': [1]})))
                    stack.enter_context(patch.object(U, 'doublet_risk', side_effect=risk))
                    model = stack.enter_context(patch.object(I, 'models', return_value=[]))
                    seeds = stack.enter_context(patch.object(I, 'seed_diagnostic'))
                    D.interactions(out / 'core', out, Mock(get=lambda url: url), SimpleNamespace(rscript='Rscript'))
                seeds.assert_not_called()
                summary = pd.read_csv(out / 'doublet_library_summary.csv')
                self.assertEqual(len(summary), len(samples))
                self.assertEqual(summary.status.iloc[0], 'not_testable_insufficient_qc_cells')
                filtered, filtered_units = model.call_args.args[3:5]
                self.assertNotIn('small', filtered.columns)
                self.assertNotIn('small', filtered_units.gsm.tolist())
                self.assertEqual(list(filtered.columns), ['good'] if include_good else [])
                self.assertEqual(list(filtered_units.columns), list(units.columns))
                self.assertTrue((out / 'doublet_seed_summary.csv').is_file())


if __name__ == '__main__':
    unittest.main()
