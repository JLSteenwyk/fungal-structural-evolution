#!/usr/bin/env python3
"""Recover original native completion evidence when its API wait handle is lost."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('launch', 'execution', 'validation', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--unavailable-original-session', type=int, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launch, execution, validation = [json.loads(p.read_text()) for p in
                                     (args.launch, args.execution, args.validation)]
    assert execution['status'] == 'exited_zero_with_receipt'
    assert execution['exit_code'] == 0 and execution['timed_out'] is False
    assert execution['receipt'] == str(args.validation) and execution['receipt_sha256'] == sha(args.validation)
    assert validation['scientific_eligibility'] is False
    assert launch['invocation_id'] == execution['invocation_id']
    for key in ('pid', 'created', 'cmdline'):
        assert launch[key] == execution['wrapper'][key]
    assert launch['actual_cgroup_limits'] == execution['actual_cgroup_limits']
    assert launch.get('original_tool_session_id') in (None, args.unavailable_original_session)
    assert sha(launch['plan']) == launch['plan_sha256']
    config_path = args.execution.with_suffix('')/'configuration.json'
    config = json.loads(config_path.read_text())
    assert all(execution[key] == value for key, value in config.items())
    native = json.loads((args.execution.with_suffix('')/'process.json').read_text())
    assert execution['child'] == native and execution['command'] == native['command']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', launch['unit'], '-o', 'json', '--no-pager'], text=True).splitlines()]
    invocation = execution['invocation_id']
    rows = [r for r in rows if invocation in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID'))]
    exact = [r for r in rows if r.get('_PID') == str(execution['wrapper']['pid'])
             and r.get('_CMDLINE') == ' '.join(execution['wrapper']['cmdline'])]
    assert len(exact) == 2
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=execution['wrapper'], invocation_id=invocation)
    terminal = {k: v for k, v in execution.items() if k not in
                ('source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope')}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    assert int(exact[0]['__REALTIME_TIMESTAMP']) < int(exact[1]['__REALTIME_TIMESTAMP'])
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == invocation and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == invocation and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    assert not any(r.get('UNIT_RESULT') not in (None, 'success')
                   or r.get('EXIT_STATUS') not in (None, '0') for r in rows)
    journal = args.execution.with_suffix('')/'recovered-original-native-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows: handle.write(json.dumps(row, sort_keys=True)+'\n')
    pins = {}
    for mapping in (execution['source_hashes'], execution['artifacts'], validation['source_hashes'], launch['source_hashes']):
        for path, digest in mapping.items(): bind(pins, path, digest)
    for path in (args.launch, args.execution, args.validation, journal, Path(launch['plan']), Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_native_receipt_and_whole_journal_without_api_wait',
                  checked_utc=datetime.now(timezone.utc).isoformat(), unit=launch['unit'], invocation_id=invocation,
                  wrapper=execution['wrapper'], child=execution['child'], native_exit_code=execution['exit_code'],
                  wrapper_initial_and_terminal_matched=True, manager_start_records=1, manager_completion_records=1,
                  unavailable_original_tool_session_id=args.unavailable_original_session,
                  original_tool_terminal_exit_code=None, original_tool_terminal_available=False,
                  validation_sha256=sha(args.validation), source_hashes=pins, scientific_eligibility=False,
                  scope='Exact original immutable launch/configuration/native identity, full original wrapper journal, '
                        'recorded native zero exit and original manager completion; all declared source/output hashes '
                        'reverified. The API wait handle is unavailable; no actual API terminal or replacement wait '
                        'is invented, and the original command is not rerun. This is native execution evidence, '
                        'not the original-tool-wait transport contract, independent matrix-statistic qualification '
                        'or biological acceptance.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'wrapper', 'child')}, indent=2))


if __name__ == '__main__': main()
