#!/usr/bin/env python3
"""Close exact original V2 GDB software controls, transport and source bindings."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    ep = Path('metadata/retained_factor_watchpoint_debugger_v2_software_execution_20261004_v1.json')
    vp = Path('metadata/retained_factor_watchpoint_debugger_v2_software_validation_20261004_v1.json')
    e = json.loads(ep.read_text()); v = json.loads(vp.read_text())
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and not e['timed_out']
    assert sha(vp) == e['receipt_sha256']
    assert v['status'] == 'passed_actual_native_watchpoint_debugger_v2_controls'
    assert v['actual_native_debugger_controls'] == 2 and v['original_factor_diagnostic_runs'] == 0
    assert v['artificial_control'] and not v['scientific_eligibility']
    assert set(v['outcomes']) == {'native_write', 'normal_exit'}
    for case, outcome in v['outcomes'].items():
        assert outcome['debugger_exit_code'] == outcome['driver_exit_code'] == 0
        assert outcome['installed_watchpoints'] == 2
        assert outcome['native_write_hit'] == (case == 'native_write')
        assert outcome['normal_exit_inspection_skipped'] == (case == 'normal_exit')
    pins = {}
    for mapping in [e['source_hashes'], e['artifacts'], v['source_hashes']]:
        for name, digest in mapping.items(): bind(pins, name, digest)
    inv = e['invocation_id']; unit = 'fungal-retained-watchpoint-v2-software-20261004-v1.service'
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [row for row in rows if inv in [row.get('_SYSTEMD_INVOCATION_ID'), row.get('USER_INVOCATION_ID')]]
    exact = [row for row in rows if row.get('_PID') == str(e['wrapper']['pid'])
        and row.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(row['_SYSTEMD_INVOCATION_ID'] == inv for row in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: val for k, val in e.items() if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and 'Started ' in row.get('MESSAGE', '')]
    ends = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and row.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows: handle.write(json.dumps(row, sort_keys=True)+'\n')
    for path in [Path(__file__), ep, vp, journal]: bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_watchpoint_debugger_v2_software_wait_zero',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_tool_session_id=54616,
        actual_tool_terminal_exit_code=0, unit=unit, invocation_id=inv, wrapper=e['wrapper'],
        original_tool_terminal_payload=dict(chunk_id='adf20b', wall_time_seconds=0.00001337,
            exit_code=0, original_token_count=44,
            output='Finished with result: success\nMain processes terminated with: code=exited/status=0\n'
                'Service runtime: 3.083s\nCPU time consumed: 3.413s\nMemory peak: 256.0K\nMemory swap peak: 0B\n'),
        entire_terminal_payload_matched=True, original_start_records=1, original_completion_records=1,
        source_hashes=pins, actual_native_debugger_controls=2, artificial_control=True,
        original_factor_diagnostic_runs=0, original_jobs_restarted=False,
        original_corruption_resolved=False, scientific_eligibility=False, biological_fits=0, gpu=False,
        scope='Original bounded actual two-case native GDB qualification; exact wrapper PID/create/command/'
              'invocation, complete initial/terminal wrapper payloads and manager start/end. Source/artifact '
              'hashes freshly verified. Native write capture and conditional normal-exit inspection both '
              'pass on artificial words. Generic wrapper kernel label is superseded by actual command. '
              'No original factor retry, memory-corruption repair, posterior or biological acceptance.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: val for k, val in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__': main()
