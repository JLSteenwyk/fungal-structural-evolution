#!/usr/bin/env python3
"""Preserve completed probes separately from their debugger's post-exit error."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    ep = Path('metadata/retained_factor_watchpoint_debugger_execution_20261004_v1.json')
    rp = Path('metadata/retained_factor_watchpoint_debugger_20261004_v1.json')
    e = json.loads(ep.read_text()); r = json.loads(rp.read_text())
    assert e['exit_code'] == 1 and e['status'] == 'failed_stage_retained' and not e['timed_out']
    assert e['receipt_sha256'] == sha(rp)
    assert r['status'] == 'retained_watchpoint_debugger_or_inferior_failure_requires_review'
    assert r['debugger_exit_code'] == 1 and not r['native_write_hit']
    assert len(r['installed_watchpoints']) == 2 and not r['artificial_control']
    plan_path = Path(r['plan']); plan = json.loads(plan_path.read_text())
    assert sha(plan_path) == r['plan_sha256']
    target_path = Path(plan['target_receipt']); target = json.loads(target_path.read_text())
    assert target['status'] == 'completed_watched_factor_timing_probe_diagnostic_v8'
    assert target['planned_groups'] == target['executed_groups'] == 40
    assert target['probe_status_counts'] == {'timed_all_declared_points_agree_only': 40}
    assert target['inputs_preserved_for_all_attempted_groups'] and not target['detected_input_mutation']
    assert not target['unattempted_group_ids'] and not target['scientific_eligibility']
    debug = Path(plan['debugger_output'])
    stdout = (debug/'gdb_stdout.log').read_text(); stderr = (debug/'gdb_stderr.log').read_text()
    assert 'exited normally]' in stdout and 'No registers.' in stderr
    assert 'Error in sourced command file:' in stderr
    assert 'PROJECT_STOP_PC=' not in stdout
    assert 'Old value = ' not in stdout and 'New value = ' not in stdout
    pins = {}
    for mapping in [e['source_hashes'], e['artifacts'], r['source_hashes'],
                    target['source_hashes'], target['artifacts']]:
        for name, digest in mapping.items(): bind(pins, name, digest)
    unit = 'fungal-retained-factor-watchpoint-debugger-20261004-v1.service'
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    inv = e['invocation_id']
    rows = [row for row in rows if inv in [row.get('_SYSTEMD_INVOCATION_ID'), row.get('USER_INVOCATION_ID')]]
    exact = [row for row in rows if row.get('_PID') == str(e['wrapper']['pid'])
        and row.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(row['_SYSTEMD_INVOCATION_ID'] == inv for row in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: v for k, v in e.items() if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and 'Started ' in row.get('MESSAGE', '')]
    ends = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and row.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows: handle.write(json.dumps(row, sort_keys=True)+'\n')
    for path in [Path(__file__), ep, rp, target_path, plan_path, journal]: bind(pins, path)
    verify(pins)
    result = dict(status='verified_completed_unprotected_probes_with_original_debugger_post_exit_failure',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_tool_session_id=91413,
        actual_tool_terminal_exit_code=1, debugger_exit_code=1, inferior_exited_normally=True,
        invocation_id=inv, unit=unit, wrapper=e['wrapper'],
        original_tool_terminal_payload=dict(chunk_id='bf0df2', wall_time_seconds=0.000016185,
            exit_code=1, original_token_count=48,
            output='Finished with result: exit-code\nMain processes terminated with: code=exited/status=1\n'
                   'Service runtime: 18min 37.665s\nCPU time consumed: 18min 39.965s\n'
                   'Memory peak: 256.0K\nMemory swap peak: 0B\n'),
        entire_wrapper_initial_and_terminal_payload_matched=True,
        original_start_records=1, original_completion_records=1,
        target_probes_checked=40, inputs_preserved=True, installed_watchpoints=2,
        native_write_observed=False, exact_native_write_instruction_identified=False,
        original_corruption_resolved=False, source_hashes=pins,
        original_jobs_restarted=False, production_tolerance_changed=False,
        biological_fits=0, scientific_eligibility=False, gpu=False,
        scope='All40actual unprotected probes completed with unchanged guards and preserved inputs. '
              'Original GDB and wrapper exit1 remain failures: sourced register inspection ran after '
              'normal inferior exit. No terminal-zero debugger qualification is claimed. Exact original '
              'wrapper/journal/source/artifact evidence preserved. Two-word watchpoint nonrecurrence '
              'under debugger does not repair prior corruption or identify its instruction. No retry.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__': main()
