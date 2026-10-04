#!/usr/bin/env python3
"""Capture the first write to two exact, unprotected source-factor words."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha, write_json
from reference_measurement_union_sources import bind, verify


def commands(metadata):
    # The inferior writes its own live addresses before raising SIGSTOP.
    python = '\n'.join(['python', 'import gdb, json',
        'watch = json.load(open(' + repr(str(metadata)) + '))',
        'assert len(watch["targets"]) == 2',
        'assert gdb.selected_inferior().pid == watch["inferior_pid"]',
        'for name, item in watch["targets"].items():',
        '    assert item["bytes"] == 8 and item["address"] % 8 == 0',
        '    gdb.execute("watch -location *(unsigned long long *)" + hex(item["address"]))',
        '    print("PROJECT_ARMED_WATCH=" + name + ":" + hex(item["address"]))', 'end'])
    return ['set pagination off', 'set confirm off', 'set debuginfod enabled off',
        'set disable-randomization off', 'set auto-load python-scripts off',
        'handle SIGSTOP stop print nopass', 'handle SIGSEGV stop print nopass',
        'run', python, 'info breakpoints', 'continue',
        'printf "\\nPROJECT_STOP_PC=%p\\n", $pc', 'thread apply all bt 40',
        'info registers', 'x/12i $pc-24', 'info sharedlibrary', 'info breakpoints']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); verify(plan['pins'])
    pins = dict(plan['pins']); bind(pins, args.plan)
    root = Path(plan['debugger_output']); root.mkdir(exist_ok=False)
    command_file = root/'gdb_commands.txt'
    with command_file.open('x') as handle:
        handle.write('\n'.join(commands(plan['watchpoint_metadata'])) + '\n')
    command = ['/usr/bin/gdb', '-nx', '--batch', '-x', str(command_file),
        '--args', *plan['target_command']]
    write_json(root/'debugger_command.json', command)
    with (root/'gdb_stdout.log').open('x') as out, (root/'gdb_stderr.log').open('x') as err:
        child = subprocess.Popen(command, stdout=out, stderr=err, env=dict(os.environ, DEBUGINFOD_URLS=''))
        proc = psutil.Process(child.pid)
        identity = dict(pid=proc.pid, created=proc.create_time(), cmdline=proc.cmdline())
        write_json(root/'debugger_process.json', identity)
        code = child.wait()
    text = (root/'gdb_stdout.log').read_text()
    armed = re.findall(r'PROJECT_ARMED_WATCH=([^:\n]+):(0x[0-9a-fA-F]+)', text)
    # Installation prints a watchpoint declaration. A hit additionally reports
    # old/new values; do not misclassify installation as a native write.
    hit = bool(re.search(r'Hardware watchpoint \d+:.*\n\s*\nOld value = .*\nNew value = ', text))
    stops = re.findall(r'PROJECT_STOP_PC=(0x[0-9a-fA-F]+)', text)
    metadata = Path(plan['watchpoint_metadata'])
    if metadata.exists():
        watch = json.loads(metadata.read_text())
        assert watch['artificial_control'] == plan['artificial_control']
        assert not watch['process_local_read_only_protection']
        assert {name: int(address, 16) for name, address in armed} == {
            name: item['address'] for name, item in watch['targets'].items()}
    target_receipt = Path(plan['target_receipt'])
    if hit and len(armed) == 2 and stops and code == 0:
        status = ('captured_native_watchpoint_positive_control' if plan['artificial_control']
            else 'captured_native_write_to_unprotected_source_factor')
    elif target_receipt.exists() and len(armed) == 2 and not hit and code == 0:
        target = json.loads(target_receipt.read_text())
        assert target['status'] == plan['target_expected_status']
        status = 'completed_unprotected_factor_watchpoint_probe_without_hit'
    else:
        status = 'retained_watchpoint_debugger_or_inferior_failure_requires_review'
    for folder in [root, Path(plan['output'])]:
        for path in folder.rglob('*'):
            if path.is_file(): bind(pins, path)
    verify(pins)
    result = dict(status=status, checked_utc=datetime.now(timezone.utc).isoformat(),
        plan=str(args.plan), plan_sha256=sha(args.plan), debugger_exit_code=code,
        debugger_identity=identity, installed_watchpoints=armed, native_write_hit=hit,
        stopped_program_counters=stops, target_receipt_present=target_receipt.exists(),
        artificial_control=plan['artificial_control'], source_hashes=pins,
        process_local_read_only_protection=False, original_jobs_restarted=False,
        production_tolerance_changed=False, scientific_eligibility=False, biological_fits=0, gpu=False,
        scope='Separate unprotected inferior with two hardware word watchpoints at the previously '
              'affected source-factor locations. First hit stops; partial evidence and native stack/PC '
              'remain. The two watched words do not cover every byte of five factors. Debugger and '
              'SIGSTOP alter execution context. Positive control is artificial, not source corruption. '
              'No automatic retries, global native-library changes, tolerance relaxation or biological acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({key: value for key, value in result.items()
        if key not in ['source_hashes', 'debugger_identity']}, indent=2))
    if status == 'retained_watchpoint_debugger_or_inferior_failure_requires_review':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
