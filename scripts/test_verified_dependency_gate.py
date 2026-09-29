#!/usr/bin/env python3
"""Exercise failed dependencies, changed pins and the process-exit/unit-state gap."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import run_after_verified_dependencies as gate


class DependencyGate(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        depplan = root/'dependency_plan.json'
        depplan.write_text('{}\n')
        launch = root/'launch.json'
        launch.write_text(json.dumps(dict(plan=str(depplan), plan_sha256=gate.sha(depplan),
            unit='test-dependency.service', pid=999999999, created=0., cmdline=['fixture'])))
        self.plan = root/'plan.json'
        self.plan.write_text(json.dumps(dict(dependencies=[str(launch)], command=['fixture-child'],
            pins={str(launch): gate.sha(launch), str(depplan): gate.sha(depplan)})))
        self.depplan = depplan
        self.success = dict(ActiveState='inactive', Result='success', ExecMainStatus='0')

    def execute(self, states):
        with patch('sys.argv', ['runner', '--plan', str(self.plan)]), \
             patch.object(gate, 'live', return_value=False), \
             patch.object(gate, 'terminal_state', side_effect=states), \
             patch.object(gate.time, 'sleep') as sleep, \
             patch.object(gate.subprocess, 'run') as child:
            gate.main()
            child.assert_called_once_with(['fixture-child'], check=True)
            return sleep.call_count

    def test_success_runs_child_once(self):
        self.assertEqual(self.execute([self.success, self.success]), 0)

    def test_failed_dependency_never_runs_child(self):
        with patch('sys.argv', ['runner', '--plan', str(self.plan)]), \
             patch.object(gate, 'live', return_value=False), \
             patch.object(gate, 'terminal_state', return_value=dict(ActiveState='failed', Result='exit-code', ExecMainStatus='3')), \
             patch.object(gate.subprocess, 'run') as child:
            with self.assertRaises(RuntimeError):
                gate.main()
            child.assert_not_called()

    def test_unit_transition_is_polled_without_running_child_early(self):
        transition = dict(ActiveState='deactivating', Result='success', ExecMainStatus='0')
        self.assertEqual(self.execute([transition, self.success, self.success]), 1)

    def test_changed_plan_never_runs_child(self):
        self.depplan.write_text('{"changed": true}\n')
        with patch('sys.argv', ['runner', '--plan', str(self.plan)]), patch.object(gate.subprocess, 'run') as child:
            with self.assertRaises(AssertionError):
                gate.main()
            child.assert_not_called()

    def test_wrong_live_command_is_rejected(self):
        import psutil
        process = psutil.Process()
        with self.assertRaises(ValueError):
            gate.live(dict(pid=process.pid, created=process.create_time(), cmdline=['different-command'], unit='fixture'))


if __name__ == '__main__':
    unittest.main()
