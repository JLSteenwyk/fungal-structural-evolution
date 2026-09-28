#!/usr/bin/env python3
"""Exercise chain recovery using disposable local processes, never real analyses."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from ancestral_chain_attempt import run_attempt, sha


class AttemptTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / 'chain'
        self.input = Path(self.temp.name) / 'input.txt'
        self.input.write_text('fixed input')

    def tearDown(self):
        self.temp.cleanup()

    def config(self, code, timeout=5):
        return dict(command=[sys.executable, '-c', code], timeout_seconds=timeout,
                    seed=123, pins={str(self.input): sha(self.input), sys.executable: sha(sys.executable)})

    def test_reuse_and_tamper_detection(self):
        config = self.config("from pathlib import Path; Path('sample.txt').write_text('sample')")
        receipt = run_attempt(self.root, config)
        self.assertEqual(run_attempt(self.root, config), receipt)
        self.assertEqual(len(list(self.root.glob('attempt-*'))), 1)
        (receipt.parent / 'sample.txt').write_text('changed')
        with self.assertRaises(AssertionError):
            run_attempt(self.root, config)

    def test_failed_attempt_is_preserved(self):
        config = self.config('raise SystemExit(2)')
        first = run_attempt(self.root, config)
        before = first.read_bytes()
        second = run_attempt(self.root, config)
        self.assertNotEqual(first, second)
        self.assertEqual(first.read_bytes(), before)

    def test_timeout_and_config_binding(self):
        config = self.config('import time; time.sleep(10)', .1)
        receipt = run_attempt(self.root, config)
        self.assertEqual(json.loads(receipt.read_text())['status'], 'timeout')
        config['seed'] = 456
        with self.assertRaises(AssertionError):
            run_attempt(self.root, config)

    def test_missing_identity_fails_closed(self):
        (self.root / 'attempt-0001').mkdir(parents=True)
        with self.assertRaises(RuntimeError):
            run_attempt(self.root, self.config('pass'))

    def test_parent_crash_live_child_and_recovery(self):
        config = self.config('import time; time.sleep(2)')
        code = ('import sys,json;sys.path.insert(0,sys.argv[1]);'
                'from ancestral_chain_attempt import run_attempt;'
                'run_attempt(sys.argv[2],json.loads(sys.argv[3]))')
        parent = subprocess.Popen([sys.executable, '-c', code,
                                   str(Path(__file__).resolve().parent), str(self.root), json.dumps(config)])
        identity = self.root / 'attempt-0001/process.json'
        child = None
        try:
            deadline = time.monotonic() + 10
            while not identity.exists():
                self.assertIsNone(parent.poll())
                self.assertLess(time.monotonic(), deadline)
                time.sleep(.02)
            child = json.loads(identity.read_text())['pid']
            parent.kill()
            parent.wait()
            with self.assertRaises(BlockingIOError):
                run_attempt(self.root, config)
            os.killpg(child, signal.SIGKILL)
            from ancestral_chain_attempt import live_group
            while live_group(child):
                self.assertLess(time.monotonic(), deadline)
                time.sleep(.02)
            recovered = run_attempt(self.root, config)
            self.assertEqual(recovered.parent.name, 'attempt-0002')
            self.assertFalse((identity.parent / 'receipt.json').exists())
            self.assertTrue(identity.exists())
        finally:
            if parent.poll() is None:
                parent.kill()
                parent.wait()
            if child is not None:
                try:
                    os.killpg(child, signal.SIGKILL)
                except ProcessLookupError:
                    pass


if __name__ == '__main__':
    unittest.main(verbosity=2)
