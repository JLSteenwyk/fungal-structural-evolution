#!/usr/bin/env python3
"""Verify an original capped watchpoint debugger wait, sources and full journal."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--unit', required=True)
    p.add_argument('--actual-tool-session', type=int, required=True)
    p.add_argument('--actual-tool-exit', type=int, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    e = json.loads(a.execution.read_text()); v = json.loads(a.receipt.read_text())
    assert e['status'] == 'exited_zero_with_receipt' and not e['timed_out']
    assert e['exit_code'] == a.actual_tool_exit == 0 and sha(a.receipt) == e['receipt_sha256']
    assert len(v['installed_watchpoints']) == 2 and v['debugger_exit_code'] == 0
    assert not v['scientific_eligibility'] and not v['original_jobs_restarted']
    if v['artificial_control']:
        assert v['status'] == 'captured_native_watchpoint_positive_control'
        assert v['native_write_hit'] and v['stopped_program_counters']
    else:
        assert v['status'] in ['captured_native_write_to_unprotected_source_factor',
            'completed_unprotected_factor_watchpoint_probe_without_hit']
    pins = {}
    for mapping in [e['source_hashes'], e['artifacts'], v['source_hashes']]:
        for path, digest in mapping.items(): bind(pins, path, digest)
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', a.unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    inv = e['invocation_id']
    rows = [row for row in rows if inv in [row.get('_SYSTEMD_INVOCATION_ID'), row.get('USER_INVOCATION_ID')]]
    exact = [row for row in rows if row.get('_PID') == str(e['wrapper']['pid'])
        and row.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(row['_SYSTEMD_INVOCATION_ID'] == inv for row in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {key: value for key, value in e.items()
        if key not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and 'Started ' in row.get('MESSAGE', '')]
    ends = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and row.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = a.execution.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows: handle.write(json.dumps(row, sort_keys=True)+'\n')
    for path in [Path(__file__), a.execution, a.receipt, journal]: bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_capped_watchpoint_debugger_wait_zero',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_tool_session_id=a.actual_tool_session,
        actual_tool_terminal_exit_code=0, invocation_id=inv, unit=a.unit, wrapper=e['wrapper'],
        entire_terminal_payload_matched=True, original_start_records=1, original_completion_records=1,
        source_hashes=pins, artificial_control=v['artificial_control'], native_write_hit=v['native_write_hit'],
        original_corruption_resolved=False, biological_fits=0, scientific_eligibility=False, gpu=False,
        scope='Actual original wait and complete initial/terminal wrapper payloads plus manager start/end '
              'verified; all retained source/artifact hashes rechecked. Positive control concerns artificial '
              'words; actual factor watchpoints cover two words and retain partial evidence at first write. '
              'No corruption repair, original-run root cause, original restart or biological inference.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__': main()
