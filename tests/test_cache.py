"""Only source downloads may cross experiment boundaries."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pipeline.cache import destination, seed_raw_cache, import_raw_sources
from pipeline.io import download, dump, sha
from pipeline.lifecycle import claim, current_source
from pipeline.sources import URLS


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.source = Path(self.temp.name) / 'old'
        self.work = Path(self.temp.name) / 'new'
        self.source.mkdir()
        self.work.mkdir()
        self.config = {'seaad_release': '2026-06-22'}

    def tearDown(self):
        self.temp.cleanup()

    def raw(self, name='mouse.rds.gz', url=None):
        path = self.source / 'raw' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'complete raw download')
        dump(str(path)+'.download.json',
             {'url': url or URLS[name], 'sha256': sha(path), 'bytes': path.stat().st_size})
        return path

    def test_raw_skips_network(self):
        self.raw()
        seed_raw_cache(self.source, self.work, self.config)
        with patch('urllib.request.urlopen', side_effect=AssertionError('no network')):
            result = download(URLS['mouse.rds.gz'], self.work/'raw/mouse.rds.gz')
        self.assertEqual(result.read_bytes(), b'complete raw download')

    def test_no_derived_or_fitted_data_imported(self):
        self.raw()
        for name in ('prepared/mouse_units.csv', 'results/panel.csv',
                     'stages/panel.json', 'raw/SEAAD/donor_pseudobulk.npz'):
            path = self.source/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('old calculation')
        self.raw('fake_result.csv', 'https://example.invalid/old-result')
        imported = seed_raw_cache(self.source, self.work, self.config)
        self.assertEqual(len(imported), 1)
        self.assertEqual({p.name for p in (self.work/'raw').iterdir()},
                         {'mouse.rds.gz', 'mouse.rds.gz.download.json'})
        self.assertFalse((self.work/'prepared').exists())

    def test_damaged_raw_stops(self):
        path = self.raw()
        path.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'damaged'):
            seed_raw_cache(self.source, self.work, self.config)

    def test_unregistered_file_not_imported(self):
        path = self.source/'raw/mouse.rds.gz'
        path.parent.mkdir()
        path.write_bytes(b'no provenance')
        self.assertEqual(seed_raw_cache(self.source, self.work, self.config), [])

    def test_different_release_not_imported(self):
        self.assertIsNone(destination('https://sea-ad-single-cell-profiling.s3.amazonaws.com/'
                                      'PFC/RNAseq/donors_objects/d_SEAAD_DFC_RNAseq_final-nuclei.2020.h5ad',
                                      '2026-06-22'))

    def test_repeat_source_list_must_match(self):
        import_raw_sources([self.source], self.work, self.config)
        with self.assertRaisesRegex(ValueError, 'cannot change'):
            import_raw_sources([], self.work, self.config)

    def test_fresh_default_rejects_old_results(self):
        (self.work/'results').mkdir()
        with self.assertRaisesRegex(ValueError, 'empty output'):
            claim(self.work, 'core')

    def test_explicit_same_run_resume(self):
        identity = claim(self.work, 'core', code='v')
        with self.assertRaisesRegex(ValueError, 'NEW output'):
            claim(self.work, 'core')
        self.assertEqual(claim(self.work, 'core', True, code='v'), identity)
        self.assertEqual(current_source(self.work, identity), identity)
        with self.assertRaisesRegex(ValueError, 'another experiment'):
            current_source(self.work, 'different')
        with self.assertRaisesRegex(ValueError, 'changed'):
            claim(self.work, 'core', True, code='changed')

    def test_resume_does_not_adopt_legacy_directory(self):
        dump(self.work/'owner.json', {'source':'official_public_downloads'})
        with self.assertRaisesRegex(ValueError, 'requires a fresh_run'):
            claim(self.work, 'core', True)

    def test_observed_layout_from_new_metadata(self):
        from pipeline import benchmark
        benchmark.load_science()
        path = self.work/'mouse_units.csv'
        path.write_text('SAMPLE,cage,genotype,sex,age\n2,6.5,5X,F,8\n1,6,WT,M,4\n')
        data = benchmark.observed_layout(path)
        self.assertEqual(data['sample'].tolist(), ['1','2'])
        self.assertEqual(data.cage.tolist(), [0,1])
        self.assertEqual(data.genotype.tolist(), [0,1])


if __name__ == '__main__':
    unittest.main()
