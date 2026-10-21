#!/usr/bin/env python3
"""Submit a prepared RepeatModeler recovery array only after verified Slurm recovery."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def terminal(states):
    allowed = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
    return bool(states) and all(state.split()[0].split('+')[0] in allowed for state in states)


def command(args):
    return subprocess.run(args, text=True, capture_output=True, check=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path)
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    if a.self_test:
        assert terminal(['COMPLETED', 'FAILED+']) and not terminal([]) and not terminal(['RUNNING'])
        print(json.dumps(dict(status='passed_repeatmodeler_scheduler_submission_literal_control'), indent=2)); return
    if a.plan is None: p.error('--plan is required unless --self-test is used')
    plan = json.loads(a.plan.read_text())
    pins = plan['pins']
    for path, expected in pins.items():
        if sha(path) != expected: raise ValueError(f'Changed pinned input: {path}')
    output = Path(plan['submission_receipt'])
    if output.exists(): raise FileExistsError('Recovery submission receipt already exists')
    while True:
        ping = command(['scontrol', 'ping'])
        if ping.returncode:
            print('waiting_for_slurm_controller', (ping.stdout + ping.stderr).strip(), flush=True); time.sleep(plan['poll_seconds']); continue
        active = command(['squeue', '--noheader', '--jobs', str(plan['original_array_job_id']), '--format=%i'])
        if active.returncode:
            print('waiting_for_original_array_query', (active.stdout + active.stderr).strip(), flush=True); time.sleep(plan['poll_seconds']); continue
        if active.stdout.strip():
            print('waiting_for_original_array_terminal', active.stdout.strip(), flush=True); time.sleep(plan['poll_seconds']); continue
        history = command(['sacct', '--noheader', '--parsable2', '--jobs', str(plan['original_array_job_id']), '--format=JobIDRaw,State'])
        rows = [line.split('|', 1)[1] for line in history.stdout.splitlines() if '|' in line and line.split('|', 1)[0].startswith(str(plan['original_array_job_id']))]
        if history.returncode or not terminal(rows):
            print('waiting_for_original_array_terminal_accounting', (history.stdout + history.stderr).strip(), flush=True); time.sleep(plan['poll_seconds']); continue
        submit = command(['sbatch', str(plan['sbatch'])])
        if submit.returncode:
            raise RuntimeError('Recovery submission failed: ' + (submit.stdout + submit.stderr).strip())
        receipt = dict(status='submitted_repeatmodeler_recovery_array_after_verified_scheduler_return', checked_utc=datetime.now(timezone.utc).isoformat(), original_array_job_id=plan['original_array_job_id'], original_array_terminal_states=rows, submission_stdout=submit.stdout.strip(), task_count=plan['task_count'], task_file=str(plan['task_file']), task_file_sha256=sha(plan['task_file']), sbatch=str(plan['sbatch']), sbatch_sha256=sha(plan['sbatch']), plan_sha256=sha(a.plan), scope='Submission custody only. Completion requires a fresh full-panel independent progress audit; this receipt does not establish repeat annotation or genome-architecture inference.')
        output.write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps(receipt, indent=2)); return

if __name__ == '__main__': main()
