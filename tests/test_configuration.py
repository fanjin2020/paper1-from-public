"""Reject incompatible releases before environment probes or result creation."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import reproduce_from_public as launcher
import run
from pipeline import benchmark, configuration, core, matched


class ConfigurationTests(unittest.TestCase):
    def write_config(self, root, **changes):
        data = configuration.load_config()
        data.update(changes)
        path = root / 'config.json'
        path.write_text(json.dumps(data), encoding='utf-8')
        return path

    def test_packaged_release_versions_match(self):
        config = configuration.load_config()
        for version in (core.VERSION, matched.VERSION, benchmark.VERSION, configuration.VERSION):
            self.assertEqual(config['analysis_version'], version)

    def test_mismatch_identifies_file_expected_and_actual_versions(self):
        with tempfile.TemporaryDirectory() as d:
            path = self.write_config(Path(d), analysis_version='obsolete-release')
            with self.assertRaises(ValueError) as error:
                configuration.load_config(path)
            self.assertIn(str(path), str(error.exception))
            self.assertIn(configuration.VERSION, str(error.exception))
            self.assertIn('obsolete-release', str(error.exception))

    def test_worker_override_does_not_write_configuration(self):
        with tempfile.TemporaryDirectory() as d:
            path = self.write_config(Path(d))
            before = path.read_bytes()
            self.assertEqual(configuration.load_config(path, workers=60)['workers'], 60)
            self.assertEqual(before, path.read_bytes())

    def test_invalid_settings_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            for change in ({'workers': 0}, {'workers': True}, {'calibration_repeats': 99},
                           {'calibration_repeats': 100.5}, {'ebayes_robust': 'false'},
                           {'dream_ebayes_robust': None}):
                with self.subTest(change=change):
                    path = self.write_config(Path(d), **change)
                    with self.assertRaises(ValueError):
                        configuration.load_config(path)

    def test_malformed_configuration_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'config.json'
            for text in ('invalid-json', '[]'):
                path.write_text(text, encoding='utf-8')
                with self.subTest(text=text), self.assertRaises(ValueError):
                    configuration.load_config(path)

    def test_launcher_checks_configuration_before_environment(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_config(root, analysis_version='obsolete-release')
            with patch.object(launcher, 'ROOT', root), \
                 patch.object(launcher, 'select_python', side_effect=AssertionError('must not probe')):
                with self.assertRaisesRegex(ValueError, 'Configuration version mismatch'):
                    launcher.main(['--check-env', '--output', str(root / 'results')])
            self.assertFalse((root / 'results').exists())

    def test_direct_runner_checks_before_commands_and_output(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_config(root, analysis_version='obsolete-release')
            with patch.object(run, 'ROOT', root), \
                 patch.object(run, 'commands', side_effect=AssertionError('must not build commands')):
                with self.assertRaisesRegex(ValueError, 'Configuration version mismatch'):
                    run.main(['--output', str(root / 'results')])
            self.assertFalse((root / 'results').exists())

    def test_core_checks_before_creating_work(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path = self.write_config(root, analysis_version='obsolete-release')
            argv = ['core', '--config', str(path), '--work', str(root / 'work'),
                    '--library', str(root / 'R_library')]
            with patch.object(core.sys, 'argv', argv), \
                 patch.object(core.platform, 'system', return_value='Linux'), \
                 contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'Configuration version mismatch'):
                    core.main()
            self.assertFalse((root / 'work').exists())
            self.assertFalse((root / 'R_library').exists())

    def test_check_env_includes_configuration(self):
        with patch.object(launcher, 'select_python', return_value=(Path('/existing/python'), {})), \
             patch.object(launcher, 'select_r', return_value=('Rscript', [], {})), \
             contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(launcher.main(['--check-env']), 0)
        self.assertIn('CONFIG ' + configuration.VERSION, output.getvalue())
        self.assertIn('"analysis_version": "' + configuration.VERSION + '"', output.getvalue())


if __name__ == '__main__':
    unittest.main()
