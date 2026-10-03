#!/usr/bin/env python3
"""Observe the exact original paired-resampling handles and native worker caps."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import psutil

from ancestral_chain_attempt import sha
from matched_predictor_branch_inputs import verify
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    path = Path('metadata/matched_predictor_resampling_plan_20261003_v1.json')
    plan = json.loads(path.read_text()); verify(plan['pins'])
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256'] == sha(path)
    handles = []; native = []; journal = []
    for index, launch in enumerate(inventory['launches']):
        record = json.loads(Path(launch).read_text()); record['launch'] = launch
        assert sha(record['plan']) == record['plan_sha256']
        process = fingerprint(record)
        if process is None:
            observed = journal_terminal(record)
        else:
            observed = live_record(record, process)
            group = next(line[3:] for line in (Path('/proc') / str(process.pid) / 'cgroup').read_text().splitlines() if line.startswith('0::'))
            folder = Path('/sys/fs/cgroup') / group.lstrip('/')
            caps = {key: (folder / key).read_text().strip() for key in ['cpu.max', 'memory.max', 'memory.swap.max']}
            cpu, memory = (8, 24) if index == 0 else (2, 16)
            assert caps == record['actual_cgroup_limits'] == {'cpu.max': str(cpu * 100000) + ' 100000', 'memory.max': str(memory * 2**30), 'memory.swap.max': '0'}
            observed['current_cgroup_limits'] = caps
            if index == 0:
                invocation = subprocess.check_output(['systemctl', '--user', 'show', record['unit'], '-p', 'InvocationID', '--value'], text=True).strip()
                assert invocation
                raw = subprocess.check_output(['journalctl', '--user', '-u', record['unit'], '-n', '30', '-o', 'json', '--no-pager'], text=True)
                for line in raw.splitlines():
                    item = json.loads(line)
                    assert item.get('USER_INVOCATION_ID', item.get('_SYSTEMD_INVOCATION_ID')) == invocation
                    journal.append({key: item[key] for key in ['MESSAGE', '__REALTIME_TIMESTAMP']})
                for child in process.children(recursive=True):
                    try:
                        created = child.create_time(); command = child.cmdline()
                        if not command or command[0] != plan['executable']: continue
                        limits = dict(address_space=child.rlimit(psutil.RLIMIT_AS), cpu=child.rlimit(psutil.RLIMIT_CPU), file=child.rlimit(psutil.RLIMIT_FSIZE))
                        assert limits['address_space'] == (2 * 2**30, 2 * 2**30)
                        assert limits['cpu'] == (300, 300) and limits['file'] == (8 * 2**20, 8 * 2**20)
                        assert command[command.index('-T') + 1] == '1' and '-te' in command and '-keep-ident' in command
                        assert child.create_time() == created
                        native.append(dict(pid=child.pid, created=created, cmdline=command, limits=limits, rss_bytes=child.memory_info().rss))
                    except (psutil.NoSuchProcess, ProcessLookupError):
                        continue
                assert len(native) <= 8
        handles.append(observed)
    completion = None
    if Path(plan['completion']).exists():
        completion = json.loads(Path(plan['completion']).read_text())
        assert sha(completion['full_hash_archive']) == completion['full_hash_archive_sha256']
    result = dict(status='verified_original_full_matched_predictor_resampling_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(path), frozen_pins_checked=len(plan['pins']),
        original_handles=handles, current_native_workers=native, recent_original_producer_journal=journal,
        expected_resampling_cases=53200, expected_native_roles=372400, accounting_closure=completion,
        scientific_eligibility=False, gpu=False, all_eight_aims_incomplete=True,
        scope='Exact original PID/create/command, current cgroup and sampled native process limits verified. Recent original producer messages are unclosed progress, not complete archive/native readback. A native worker exiting during observation is omitted, never treated as failed or restarted. No process/source changes or GPU use.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['original_handles', 'current_native_workers', 'recent_original_producer_journal', 'accounting_closure']}))
    print('current_native_workers', len(native))
    for row in journal:
        if row['MESSAGE'].startswith('matched_predictor_resample'): print(row['MESSAGE'])


if __name__ == '__main__':
    main()
