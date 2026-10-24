"""Version, code identity and coverage regression checks without network access."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run
import reproduce_from_public as launcher
from pipeline import core, coverage, environment, lifecycle


class EnvironmentTests(unittest.TestCase):
    def fixture(self, root):
        for name in ('pipeline', 'analysis', 'tests', 'inputs'):
            (root / name).mkdir()
        for name in ('run.py', 'config.json'):
            (root / name).write_text('initial', encoding='utf-8')
        (root / 'requirements.txt').write_text('numpy==2.2.6\n', encoding='utf-8')
        (root / 'inputs/environment.json').write_text('{}', encoding='utf-8')

    def record(self):
        return {'python': {'python': '3.12.2', 'packages': {'numpy': '2.2.6'}, 'executable': 'a'},
                'R': {'R': '4.5.3', 'packages': {'limma': '3.66.0'}, 'paths': {'limma': 'a'}}}

    def test_environment_file_changes_code_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.fixture(root)
            identity = lifecycle.code_identity(root)
            (root / 'inputs/environment.json').write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError, 'changed during execution'):
                lifecycle.assert_code(root, identity)

    def test_tests_are_part_of_release_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.fixture(root)
            identity = lifecycle.code_identity(root)
            (root / 'tests/test_new.py').write_text('pass')
            self.assertNotEqual(identity, lifecycle.code_identity(root))

    def test_launcher_is_part_of_release_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.fixture(root)
            identity = lifecycle.code_identity(root)
            (root / 'reproduce_from_public.py').write_text('pass')
            self.assertNotEqual(identity, lifecycle.code_identity(root))

    def test_version_drift_rejected_without_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'software.json'
            record = self.record()
            environment.freeze_environment(path, record)
            before = path.read_bytes()
            record['python']['packages']['numpy'] = '2.5.3'
            with self.assertRaisesRegex(ValueError, 'Software versions changed'):
                environment.freeze_environment(path, record)
            self.assertEqual(before, path.read_bytes())

    def test_r_version_drift_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'software.json'
            record = self.record()
            environment.freeze_environment(path, record)
            record['R']['packages']['limma'] = 'other'
            with self.assertRaisesRegex(ValueError, 'Software versions changed'):
                environment.freeze_environment(path, record)

    def test_installation_location_is_not_science(self):
        a, b = self.record(), self.record()
        b['python']['executable'] = 'b'
        b['R']['paths']['limma'] = 'b'
        self.assertEqual(environment.scientific_identity(a), environment.scientific_identity(b))

    def test_transitive_r_version_is_part_of_identity(self):
        a, b = self.record(), self.record()
        a['R']['loaded_namespaces'] = {'Rcpp': '1.1.2'}
        b['R']['loaded_namespaces'] = {'Rcpp': 'other'}
        self.assertNotEqual(environment.scientific_identity(a), environment.scientific_identity(b))

    def test_exact_python_pins(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.fixture(root)
            with patch.object(environment.sys, 'version_info', (3,12,0)), patch.object(environment.importlib, 'import_module', return_value=SimpleNamespace(__version__='2.5.3')):
                with self.assertRaisesRegex(ValueError, 'expected 2.2.6, found 2.5.3'):
                    environment.python_snapshot(root)

    def test_r_probe_checks_failure(self):
        with patch.object(environment.subprocess, 'run', return_value=SimpleNamespace(returncode=1, stdout='', stderr='missing limma')):
            with self.assertRaisesRegex(ValueError, 'missing limma'):
                environment.r_snapshot('Rscript', [])

    def test_r_probe_parses_marked_record(self):
        result = SimpleNamespace(returncode=0, stderr='', stdout='message\nPAPER1_ENVIRONMENT\t{"R":"4.5.3"}\n')
        with patch.object(environment.subprocess, 'run', return_value=result):
            self.assertEqual(environment.r_snapshot('Rscript', []), {'R': '4.5.3'})

    def test_full_coverage_gate_before_subprocess(self):
        from pipeline import coverage
        missing={'complete_paper':False,'missing':[{'id':'example','description':'missing implementation'}]}
        with patch.object(coverage,'report',return_value=missing), \
             patch.object(run.subprocess, 'run', side_effect=AssertionError('must not run')):
            with self.assertRaisesRegex(ValueError, 'FULL-PAPER REPRODUCTION NOT READY'):
                run.main(['--require-complete-paper'])

    def test_default_launcher_cannot_silently_run_partial(self):
        from pipeline import coverage
        missing={'complete_paper':False,'missing':[{'id':'example','description':'missing implementation'}]}
        with patch.object(coverage,'report',return_value=missing), \
             patch.object(launcher, 'select_python', side_effect=AssertionError('must not probe')):
            with self.assertRaisesRegex(ValueError, 'FULL-PAPER REPRODUCTION NOT READY'):
                launcher.main([])

    def test_coverage_is_read_only(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'not_created'
            with contextlib.redirect_stdout(io.StringIO()) as text:
                self.assertEqual(run.main(['--coverage', '--output', str(out)]), 0)
            self.assertTrue(json.loads(text.getvalue())['complete_paper'])
            self.assertFalse(out.exists())

    def test_primary_only_is_explicit(self):
        self.assertFalse(launcher.parser().parse_args([]).primary_only)
        self.assertTrue(launcher.parser().parse_args(['--primary-only']).primary_only)

    def test_bad_python_candidate_never_selected(self):
        with patch.object(launcher, 'python_candidates', return_value=[Path('/bad/python')]), \
             patch.object(launcher.subprocess, 'run', return_value=SimpleNamespace(returncode=1, stdout='', stderr='bad')):
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(RuntimeError, 'No compatible'):
                launcher.select_python()

    def test_good_python_candidate(self):
        result = SimpleNamespace(returncode=0, stderr='', stdout='PAPER1_PYTHON\t{"python":"3.12.2"}\n')
        with patch.object(launcher, 'python_candidates', return_value=[Path('/good/python')]), \
             patch.object(launcher.subprocess, 'run', return_value=result):
            path, record = launcher.select_python()
            self.assertEqual(record['python'], '3.12.2')
            self.assertEqual(path, Path('/good/python'))

    def test_explicit_python_path_is_not_resolved(self):
        with tempfile.TemporaryDirectory() as d:
            exe = Path(d) / 'venv' / 'bin' / 'python'
            exe.parent.mkdir(parents=True)
            exe.touch()
            with patch.object(Path, 'resolve', side_effect=AssertionError('must retain venv path')):
                self.assertEqual(launcher.python_candidates(str(exe)), [exe.absolute()])

    def test_nested_venv_is_discovered(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / 'paper1-from-publicV1'
            exe = root / 'paper1-from-public' / '.raw-reproduction-venv' / 'bin' / 'python'
            exe.parent.mkdir(parents=True)
            exe.touch()
            with patch.object(launcher, 'ROOT', root), \
                 patch.object(Path, 'resolve', side_effect=AssertionError('must retain venv path')):
                candidates = launcher.python_candidates()
            self.assertIn(exe.absolute(), candidates)
            self.assertEqual(candidates.count(exe.absolute()), 1)

    def test_venvs_sharing_python_target_remain_distinct(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / 'paper1-from-publicV1'
            target = Path(sys.executable).resolve()
            executables = [root / name / 'bin' / 'python'
                           for name in ('.raw-reproduction-venv', '.venv')]
            for exe in executables:
                exe.parent.mkdir(parents=True)
                try:
                    exe.symlink_to(target)
                except (OSError, NotImplementedError) as exc:
                    self.skipTest('OS does not permit executable symlinks: ' + str(exc))
            with patch.object(launcher, 'ROOT', root), patch.object(launcher.sys, 'executable', str(target)):
                candidates = launcher.python_candidates()
                self.assertIn(target, candidates)
                for exe in executables:
                    self.assertIn(exe, candidates)
                    self.assertEqual(launcher.python_candidates(str(exe)), [exe])
                bad = SimpleNamespace(returncode=1, stdout='', stderr='base version mismatch')
                good = SimpleNamespace(returncode=0, stderr='', stdout='PAPER1_PYTHON\t{"python":"3.12.2"}\n')
                with patch.object(launcher.subprocess, 'run', side_effect=[bad, good]) as probe, \
                     contextlib.redirect_stdout(io.StringIO()):
                    selected, _ = launcher.select_python()
                self.assertEqual(selected, executables[0])
                self.assertEqual(probe.call_args_list[1].args[0][0], str(executables[0]))

    def test_related_software_environment_is_discovered(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / 'paper1-from-publicV1'
            root.mkdir()
            exe = Path(d) / 'paper1_public_evidence' / 'paper3_mechanism_source_example' / '.venv' / 'bin' / 'python'
            exe.parent.mkdir(parents=True)
            exe.touch()
            with patch.object(launcher, 'ROOT', root):
                self.assertIn(exe, launcher.python_candidates())

    def test_broken_executable_does_not_hide_compatible_python(self):
        good = SimpleNamespace(returncode=0, stderr='', stdout='PAPER1_PYTHON\t{"python":"3.12.2"}\n')
        with patch.object(launcher, 'python_candidates', return_value=[Path('/bad/python'), Path('/good/python')]), \
             patch.object(launcher.subprocess, 'run', side_effect=[OSError('cannot execute'), good]), \
             contextlib.redirect_stdout(io.StringIO()):
            path, _ = launcher.select_python()
            self.assertEqual(path, Path('/good/python'))

    def test_timeout_does_not_hide_compatible_python(self):
        good = SimpleNamespace(returncode=0, stderr='', stdout='PAPER1_PYTHON\t{"python":"3.12.2"}\n')
        error = launcher.subprocess.TimeoutExpired('python', 180)
        with patch.object(launcher, 'python_candidates', return_value=[Path('/bad/python'), Path('/good/python')]), \
             patch.object(launcher.subprocess, 'run', side_effect=[error, good]), \
             contextlib.redirect_stdout(io.StringIO()):
            path, _ = launcher.select_python()
            self.assertEqual(path, Path('/good/python'))

    def test_corrupt_probe_record_is_rejected(self):
        result = SimpleNamespace(returncode=0, stderr='', stdout='PAPER1_PYTHON\tINVALID\n')
        with patch.object(launcher, 'python_candidates', return_value=[Path('/bad/python')]), \
             patch.object(launcher.subprocess, 'run', return_value=result), \
             contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError, 'No compatible'):
                launcher.select_python()

    def stage_fixture(self, base):
        root, work = base / 'source', base / 'work'
        root.mkdir()
        work.mkdir()
        self.fixture(root)
        (work / 'stages').mkdir()
        (work / 'fresh_run.json').write_text(json.dumps({'code': lifecycle.code_identity(root)}))
        (work / 'owner.json').write_text(json.dumps({'signature': 'fixture'}))
        return root, work

    def test_stage_rejects_empty_receipt(self):
        with tempfile.TemporaryDirectory() as d:
            root, work = self.stage_fixture(Path(d))
            (work / 'stages/example.json').write_text(json.dumps({'signature': 'fixture', 'outputs': {}}))
            with patch.object(core, 'ROOT', root), self.assertRaises(ValueError):
                core.stage('example', lambda: self.fail('must not run'), ['result.csv'], work, True)

    def test_stage_rejects_receipt_path_escape(self):
        with tempfile.TemporaryDirectory() as d:
            root, work = self.stage_fixture(Path(d))
            receipt = {'signature': 'fixture', 'outputs': {'../outside.csv': 'invalid'}}
            (work / 'stages/example.json').write_text(json.dumps(receipt))
            with patch.object(core, 'ROOT', root), self.assertRaisesRegex(ValueError, 'outside'):
                core.stage('example', lambda: self.fail('must not run'), ['result.csv'], work, True)

    def test_stage_checks_source_before_work(self):
        with tempfile.TemporaryDirectory() as d:
            root, work = self.stage_fixture(Path(d))
            (root / 'run.py').write_text('changed')
            with patch.object(core, 'ROOT', root), self.assertRaisesRegex(ValueError, 'changed during execution'):
                core.stage('example', lambda: self.fail('must not run'), ['result.csv'], work, False)

    def test_stage_never_receipts_work_after_source_mutation(self):
        with tempfile.TemporaryDirectory() as d:
            root, work = self.stage_fixture(Path(d))
            def compute():
                (work / 'result.csv').write_text('result')
                (root / 'run.py').write_text('changed')
            with patch.object(core, 'ROOT', root), contextlib.redirect_stdout(io.StringIO()), \
                 self.assertRaisesRegex(ValueError, 'changed during execution'):
                core.stage('example', compute, ['result.csv'], work, False)
            self.assertFalse((work / 'stages/example.json').exists())


if __name__ == '__main__':
    unittest.main()
