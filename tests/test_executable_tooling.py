"""Compile-count, cache invalidation, and failure handling for suite executables."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import test_executable as runner
from build_state import LibraryBuild


class SuiteExecutableTests(unittest.TestCase):
    def setUp(self):
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory())).resolve()
        self.programs = []
        for folder in ['alpha', 'beta']:
            path = self.root / 'tests' / folder / 'test_same.mojo'
            path.parent.mkdir(parents=True)
            path.write_text('def main() raises:\n    print("ok")\n')
            self.programs.append(path)
        self.snapshot = LibraryBuild(self.root / 'build', dict(build_id='one', compiler='mojo', artifacts={}))
        self.calls = []
        self.stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        self.stack.enter_context(patch.object(subprocess, 'run', side_effect=self.command))

    def command(self, command, **kwargs):
        self.calls.append(command)
        if command[:2] == ['mojo', 'build']:
            Path(command[command.index('-o') + 1]).write_bytes(b'one executable')
        return subprocess.CompletedProcess(command, 0, 'ok\n', '')

    def prepare(self, **kwargs):
        return runner.prepare(self.programs, self.snapshot, self.root, **kwargs)

    def test_all_cases_compile_once_and_reruns_do_not_compile(self):
        executable = self.prepare()
        self.assertFalse(executable.reused)
        for p in self.programs:
            self.assertEqual(executable.execute(p, self.root)['returncode'], 0)
        cached = self.prepare()
        cached.execute(self.programs[1], self.root)
        self.assertTrue(cached.reused)
        self.assertEqual(sum(c[:2] == ['mojo', 'build'] for c in self.calls), 1)
        self.assertEqual(len({c[0] for c in self.calls if c[0] != 'mojo'}), 1)
        self.assertEqual(self.calls[-1][-1], '1')
        self.assertNotIn(str(self.root / 'src'), self.calls[0])

    def test_same_basename_has_distinct_modules(self):
        rendered = runner.render(self.programs, self.root)
        self.assertIn('_test_cases.tests.alpha.test_same import _suite_entry_point as _case_0', rendered)
        self.assertIn('_test_cases.tests.beta.test_same import _suite_entry_point as _case_1', rendered)

    def test_entry_point_adapter_preserves_strings_and_other_definitions(self):
        source = '\"\"\"def main(): example only\"\"\"\n# def main(): comment\ndef main() raises:\n    print("main")\n'
        adapted = runner.module_source(source)
        self.assertIn('def _suite_entry_point() raises:', adapted)
        self.assertIn('def main(): example only', adapted)
        self.assertIn('# def main(): comment', adapted)
        self.assertIn('print("main")', adapted)

    def test_test_helper_library_and_flags_each_invalidate_cache(self):
        self.prepare()
        self.programs[0].write_text('def main() raises:\n    print("changed")\n')
        self.assertFalse(self.prepare().reused)
        support = self.root / 'tests/support/helper.mojo'
        support.parent.mkdir(); support.write_text('# helper\n')
        self.assertFalse(self.prepare().reused)
        self.snapshot.manifest['build_id'] = 'two'
        self.assertFalse(self.prepare().reused)
        self.assertFalse(self.prepare(flags=('-O0',)).reused)
        self.assertEqual(sum(c[:2] == ['mojo', 'build'] for c in self.calls), 5)

    def test_replaced_binary_is_never_reused(self):
        executable = self.prepare()
        executable.binary.write_bytes(b'replaced')
        self.assertFalse(self.prepare().reused)

    def test_failure_has_no_successful_manifest_and_does_not_run_tests(self):
        with patch.object(subprocess, 'run', return_value=subprocess.CompletedProcess([], 1, '', 'compile failure')):
            with self.assertRaises(runner.CompilationFailed):
                self.prepare()
        self.assertEqual(list((self.root / '.cache').rglob('manifest.json')), [])
        self.assertEqual(list((self.root / '.cache').rglob('runtime')), [])

    def test_source_change_during_compilation_is_rejected(self):
        original = self.command
        def changed(command, **kwargs):
            result = original(command, **kwargs)
            self.programs[0].write_text('# modified during compilation\n')
            return result
        with patch.object(subprocess, 'run', side_effect=changed):
            with self.assertRaisesRegex(RuntimeError, 'changed during compilation'):
                self.prepare()
        self.assertEqual(list((self.root / '.cache').rglob('manifest.json')), [])

    def test_private_imports_only_for_kernel_configuration(self):
        library = self.prepare()
        kernel = self.prepare(kernel=True)
        self.assertNotIn(str(self.snapshot.path / 'kernel'), library.record['command'])
        self.assertIn(str(self.snapshot.path / 'kernel'), kernel.record['command'])
        self.assertLess(kernel.record['command'].index(str(self.snapshot.path)), kernel.record['command'].index(str(self.snapshot.path / 'kernel')))

    def test_compiler_threads_are_bounded_and_part_of_the_cache(self):
        executable = self.prepare()
        command = executable.record['command']
        self.assertEqual(command[command.index('--num-threads') + 1], '2')
        self.assertFalse(self.prepare(compiler_threads=1).reused)
        with self.assertRaisesRegex(ValueError, 'must be positive'):
            self.prepare(compiler_threads=0)

    def test_killed_compiler_reports_its_exit_status(self):
        with patch.object(subprocess, 'run', return_value=subprocess.CompletedProcess([], -9, '', '')):
            with self.assertRaisesRegex(runner.CompilationFailed, 'status -9'):
                self.prepare()

    def test_timeout_is_a_failure(self):
        with patch.object(subprocess, 'run', side_effect=subprocess.TimeoutExpired(['test'], 1)):
            self.assertIsNone(runner.run(['test'], self.root, 1)['returncode'])


if __name__ == '__main__':
    unittest.main()
