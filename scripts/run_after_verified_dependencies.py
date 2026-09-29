#!/usr/bin/env python3
"""Run a pinned command after exact dependencies terminate successfully."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import psutil
from screen_duplication_alignment_reuse import sha


def live(launch):
    try:
        process = psutil.Process(launch['pid'])
        if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
            return False
        if process.cmdline() != launch['cmdline']:
            raise ValueError('Dependency command changed: '+launch['unit'])
        return True
    except psutil.NoSuchProcess:
        return False


def terminal_state(unit):
    return dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    dependencies = []
    for path in plan['dependencies']:
        assert str(path) in bindings
        launch = json.loads(Path(path).read_text())
        assert sha(launch['plan']) == launch['plan_sha256']
        if launch['plan'] in bindings:
            assert bindings[launch['plan']] == launch['plan_sha256']
        bindings[launch['plan']] = launch['plan_sha256']
        dependencies.append(launch)
    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pinned dependency: '+path)
    verify()
    remaining = dependencies.copy()
    while remaining:
        pending = []
        for launch in remaining:
            if live(launch):
                pending.append(launch)
                continue
            state = terminal_state(launch['unit'])
            if state['ActiveState'] in ['active', 'activating', 'deactivating', 'reloading']:
                pending.append(launch)
                continue
            if state != dict(ActiveState='inactive', Result='success', ExecMainStatus='0'):
                raise RuntimeError('Unsuccessful dependency: '+launch['unit']+' '+str(state))
        remaining = pending
        if remaining:
            print('waiting_for_successful_dependencies', [l['unit'] for l in remaining], flush=True)
            time.sleep(30)
    verify()
    # Recheck every terminal state immediately before the child runs. The wrapper
    # keeps its own PID/command unchanged, including while its child is active.
    for launch in dependencies:
        assert not live(launch)
        assert terminal_state(launch['unit']) == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    command = plan['command']
    assert isinstance(command, list) and command and all(isinstance(v, str) for v in command)
    print('running_verified_dependent_command', command, flush=True)
    subprocess.run(command, check=True)
    verify()


if __name__ == '__main__':
    main()
