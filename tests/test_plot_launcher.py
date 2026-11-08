"""Checks for the existing-environment publication figure launcher."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('plot_launcher', ROOT / 'paper1_plot_figures.py')
L = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(L)


class PlotLauncherTests(unittest.TestCase):
    def test_new_short_output_beside_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'repository'
            root.mkdir()
            first = L.create_output(root, None)
            second = L.create_output(root, None)
            self.assertEqual(first.parent, root.parent)
            self.assertTrue(first.name.startswith('p1_fig_'))
            self.assertLess(len(first.name), 40)
            self.assertNotEqual(first, second)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            sentinel = path / 'keep.txt'
            sentinel.write_text('unchanged')
            with self.assertRaisesRegex(RuntimeError, 'already exists'):
                L.create_output(path, path)
            self.assertEqual(sentinel.read_text(), 'unchanged')

    def test_explicit_python_does_not_fall_back(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(L.sys, 'executable', '/irrelevant/python'):
                self.assertEqual(list(L.candidates(root, root / 'missing')), [])

    def test_missing_bundle_stops_before_environment_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(L, 'select_python') as select:
                with self.assertRaisesRegex(RuntimeError, 'complete publication_figures'):
                    L.main(['--project-root', directory])
                select.assert_not_called()

    def test_check_env_does_not_create_output(self):
        snapshot = dict(packages={}, font='Arial')
        with patch.object(L, 'select_python', return_value=(Path('python'), snapshot)):
            with patch.object(L, 'create_output') as create:
                self.assertEqual(L.main(['--project-root', str(ROOT), '--check-env']), 0)
                create.assert_not_called()


if __name__ == '__main__':
    unittest.main()
