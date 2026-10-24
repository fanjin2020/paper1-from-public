"""Supplementary identities, raw-only imports and restart safeguards."""
import builtins
import json
from pathlib import Path
import symtable
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from pipeline.io import sha
from pipeline.supplementary import common as C, driver as D, layouts as L
from pipeline.supplementary import interactions as I, interaction_math as U
from pipeline.supplementary.conditional_simulation import simulate_designs


class SupplementaryTests(unittest.TestCase):
    def test_simulation_order_record_preserves_order_and_detects_changes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); core = root/'core'; out = root/'simulation'
            (core/'prepared').mkdir(parents=True); out.mkdir()
            data = pd.DataFrame({'SAMPLE': [652, 211, 351], 'cage': ['6.5', '2.1', '3.5'],
                                 'genotype': [0, 1, 0], 'age': [6, 12, 6],
                                 'male': [1, 0, 1], 'slide': ['s2', 's1', 's2']})
            data.to_csv(core/'prepared/mouse_units.csv', index=False)
            original = data.copy(deep=True)
            record = D.record_simulation_order(core, out, data)
            pd.testing.assert_frame_equal(data, original)
            saved = pd.read_csv(out/'simulation_animal_order.csv')
            self.assertEqual(saved.SAMPLE.tolist(), [652, 211, 351])
            self.assertEqual(saved.draw_position.tolist(), [1, 2, 3])
            self.assertEqual(D.record_simulation_order(core, out, data), record)
            with self.assertRaisesRegex(ValueError, 'input/order changed'):
                D.record_simulation_order(core, out, data.iloc[::-1])
            (out/'simulation_animal_order.csv').write_text('damaged')
            with self.assertRaisesRegex(ValueError, 'missing or damaged'):
                D.record_simulation_order(core, out, data)

    def test_simulation_order_rejects_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            data = pd.DataFrame({'SAMPLE': [211, 211], 'cage': ['2.1', '2.1'],
                                 'genotype': [0, 1], 'age': [6, 6],
                                 'male': [0, 1], 'slide': ['s1', 's1']})
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                D.record_simulation_order(root, root, data)

    def test_no_unresolved_module_globals(self):
        allowed = set(dir(builtins)) | {'__name__', '__file__', '__package__'}
        def descend(table):
            return [table] + [sub for child in table.get_children() for sub in descend(child)]
        for file in Path(D.__file__).parent.glob('*.py'):
            root = symtable.symtable(file.read_text(encoding='utf-8'), str(file), 'exec')
            missing = {s.get_name() for table in descend(root)[1:] for s in table.get_symbols()
                       if s.is_referenced() and s.is_global() and s.get_name() not in allowed
                       and s.get_name() not in root.get_identifiers()}
            self.assertFalse(missing, f'{file}: {missing}')

    def test_raw_cache_url_and_sha(self):
        url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/matrix/GSE1297_series_matrix.txt.gz'
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); cache = root / 'cache'; cache.mkdir()
            f = cache / 'download.gz'; f.write_bytes(b'public bytes')
            f.with_name(f.name + '.download.json').write_text(json.dumps({'url': url, 'sha256': sha(f)}))
            dl = C.Download(root / 'new', [cache], offline=True)
            with patch.object(C, 'download', wraps=C.download) as download:
                result = dl.get(url)
                self.assertEqual(result.read_bytes(), f.read_bytes())
                download.assert_called_once()
            with self.assertRaisesRegex(C.Blocked, 'Offline'):
                dl.get(url.replace('1297', '48350'))

    def test_unverified_matrix_is_not_imported(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); cache = root / 'cache'; cache.mkdir()
            (cache / 'GSE1297_series_matrix.txt.gz').write_bytes(b'no receipt')
            dl = C.Download(root / 'new', [cache], offline=True)
            with self.assertRaisesRegex(C.Blocked, 'Offline'):
                dl.get('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/matrix/GSE1297_series_matrix.txt.gz')

    def test_planned_bh_preserves_missing_slots(self):
        p = np.array([.001, np.nan, .025, .9])
        q = C.reserved_bh(p)
        expected = false_discovery_control(np.nan_to_num(p, nan=1.))
        self.assertTrue(np.isnan(q[1]))
        np.testing.assert_allclose(q[[0, 2, 3]], expected[[0, 2, 3]])
        with self.assertRaises(C.Blocked):
            C.reserved_bh([float('inf')])

    def test_empty_interaction_retains_all_endpoint_slots(self):
        endpoints = pd.DataFrame({'gene_symbol': ['A', 'B', 'C']})
        grid = I.fixed_grid({'endpoints': endpoints, 'frames': []})
        self.assertEqual(len(grid), 2 * 6 * 3)
        self.assertTrue(grid.P_Value.isna().all())
        self.assertTrue(grid.within_variant_family_size.eq(18).all())
        self.assertTrue(grid.joint_family_size.eq(36).all())

    def test_six_r_templates_use_exact_comparator_and_second_rho(self):
        for comparator in I.COMPARATORS:
            script, arms = I.common_background_engine(comparator)
            self.assertIn('rho<-rho_second', script)
            self.assertIn('keep<-ka&ko&koa&koo', script)
            self.assertIn('status=42', script)
            self.assertEqual(arms, ['Ast_minus_' + comparator + '_interaction_second_rho'])

    def test_stage_partial_not_reused_completed_is_verified(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            partial = root / 'cognition/attempt_001'; partial.mkdir(parents=True)
            (partial / 'result.csv').write_text('unfinished')
            def build(out):
                (out / 'result.csv').write_text('fresh')
            final = D.run_stage('cognition', root, build)
            self.assertNotEqual(final, partial)
            again = D.run_stage('cognition', root, lambda _: self.fail('must use verified same-run receipt'))
            self.assertEqual(final, again)
            (final / 'result.csv').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'damaged'):
                D.run_stage('cognition', root, build)

    def test_current_files_cannot_escape(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); (root / 's').mkdir(); (root / 'old.csv').write_text('old')
            with self.assertRaises(C.Blocked):
                C.CurrentFiles(root / 's').get('../old.csv')

    def test_designs_use_supplied_animals_not_embedded_layout(self):
        u = L.synthetic_units(4, .5, .5, 8)
        u['age'] = np.tile([4, 6, 8, 12], 18)
        rng = np.random.default_rng(124)
        u['male'] = rng.permutation(np.tile([0, 1, 1, 0], 18))
        u['slide'] = rng.permutation(np.tile([0, 1, 2], 24))
        layouts, registry = L.build_layouts(3252026, u)
        self.assertEqual(registry.loc[0, 'n_animals'], 72)
        self.assertEqual(registry.loc[0, 'n_cages'], 18)
        pd.testing.assert_series_equal(layouts['actual_linear'][0].cage, u.cage)

    def test_conditional_simulation_deterministic(self):
        rng = np.random.default_rng(23)
        g = np.tile([0., 1.], 12); c = np.repeat(np.arange(6), 4)
        base = np.column_stack([np.ones(24), g, rng.normal(size=24)])
        fixed = np.column_stack([base, pd.get_dummies(c, drop_first=True).to_numpy(float)])
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            simulate_designs({'example': (base, fixed, c)}, out, 20, 17)
            a = pd.read_csv(out / '322_design_simulation.csv')
            simulate_designs({'example': (base, fixed, c)}, out, 20, 17)
            pd.testing.assert_frame_equal(a, pd.read_csv(out / '322_design_simulation.csv'))
            self.assertEqual(len(a), 24)
            self.assertTrue(a.valid.eq(20).all())


if __name__ == '__main__':
    unittest.main()
