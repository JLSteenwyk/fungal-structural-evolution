#!/usr/bin/env python3
"""Crash only an isolated test controller and verify real BAli-Phy recovery."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil
from ancestral_chain_attempt import run_attempt, sha, live_group, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    for name, digest in plan['pins'].items():
        assert sha(name) == digest, name
    out = Path(plan['output']).resolve()
    out.mkdir(parents=True, exist_ok=False)
    root = out / 'chain'
    config = plan['config']
    code = ('import sys,json;sys.path.insert(0,sys.argv[1]);'
            'from ancestral_chain_attempt import run_attempt;'
            'run_attempt(sys.argv[2],json.loads(sys.argv[3]))')
    controller = subprocess.Popen([sys.executable, '-c', code,
        str(Path(__file__).resolve().parent), str(root), json.dumps(config)])
    identity_path = root / 'attempt-0001/process.json'
    identity = None
    try:
        deadline = time.monotonic() + 120
        while True:
            assert controller.poll() is None, 'Controller ended before crash check'
            assert time.monotonic() < deadline, 'No active sampler log before deadline'
            if identity_path.exists():
                identity = json.loads(identity_path.read_text())
                logs = list(identity_path.parent.glob('recovery-check-*/C1.log'))
                if logs and logs[0].stat().st_size > 0:
                    break
            time.sleep(.1)
        child = psutil.Process(identity['pid'])
        assert child.create_time() == identity['created']
        assert child.cmdline() == identity['command']
        assert os.getpgid(child.pid) == identity['pgid']
        lock = root / 'attempt.lock'
        retained_fds = []
        for fd in Path('/proc/%s/fd' % child.pid).iterdir():
            try:
                if os.path.samefile(fd, lock):
                    retained_fds.append(fd.name)
            except FileNotFoundError:
                pass
        assert retained_fds, 'Sampler closed inherited lock descriptor'
        controller.kill()
        controller.wait()
        blocked = False
        try:
            run_attempt(root, config)
        except BlockingIOError:
            blocked = True
        assert blocked, 'Duplicate attempt was not blocked'
        assert child.create_time() == identity['created'] and child.is_running()
        os.killpg(child.pid, signal.SIGKILL)
        while live_group(child.pid):
            assert time.monotonic() < deadline
            time.sleep(.1)
        interrupted = {str(p.relative_to(identity_path.parent)): sha(p)
                       for p in identity_path.parent.rglob('*') if p.is_file()}
        result_path = run_attempt(root, config)
        assert result_path.parent.name == 'attempt-0002'
        result = json.loads(result_path.read_text())
        assert result['exit_code'] == 0 and result['status'] == 'exited_zero_pending_scientific_validation'
        assert run_attempt(root, config) == result_path
        for name, digest in interrupted.items():
            assert sha(identity_path.parent / name) == digest
        assert not (identity_path.parent / 'receipt.json').exists()
        logs = list(result_path.parent.glob('recovery-check-*/C1.log'))
        assert len(logs) == 1
        rows = list(csv.DictReader(logs[0].open(), delimiter='\t'))
        assert [int(r['iter']) for r in rows] == list(range(21))
        for row in rows:
            prior, likelihood, posterior = [float(row[k]) for k in ['prior', 'likelihood', 'posterior']]
            assert all(math.isfinite(x) for x in [prior, likelihood, posterior])
            assert abs(prior + likelihood - posterior) < 1e-7
        receipt = dict(status='passed_isolated_baliphy_attempt_recovery',
            plan_sha256=sha(args.plan), retained_lock_fds=retained_fds,
            duplicate_launch_blocked=blocked, interrupted_artifacts=interrupted,
            recovered_attempt_receipt=str(result_path), recovered_receipt_sha256=sha(result_path),
            iterations_checked=21,
            scope='Installed BAli-Phy, one 77-tip input, controller crash and fresh whole-chain rerun. '
                  'Not native checkpoint continuation, representative resource estimate, '
                  'posterior convergence, or scientific qualification of ancestral samples.')
        write_json(out / 'receipt.json', receipt)
        print(json.dumps(receipt), flush=True)
    finally:
        if controller.poll() is None:
            controller.kill()
            controller.wait()
        if identity is not None and live_group(identity['pgid']):
            try:
                process = psutil.Process(identity['pid'])
                if process.create_time() == identity['created']:
                    os.killpg(identity['pgid'], signal.SIGKILL)
            except (ProcessLookupError, psutil.NoSuchProcess):
                pass


if __name__ == '__main__':
    main()
