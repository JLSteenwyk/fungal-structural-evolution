#!/usr/bin/env python3
"""Close actual preparer/debugger executions using their original journal payloads."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['preparation', 'protected-debugger'], required=True)
    parser.add_argument('--execution', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--unit', required=True)
    parser.add_argument('--original-tool-session', type=int, required=True)
    parser.add_argument('--original-tool-exit', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    execution = json.loads(args.execution.read_text())
    receipt = json.loads(args.receipt.read_text())
    assert execution['exit_code'] == args.original_tool_exit == 0
    assert execution['status'] == 'exited_zero_with_receipt' and not execution['timed_out']
    assert sha(args.receipt) == execution['receipt_sha256']
    pins = {}
    for mapping in [execution['source_hashes'], execution['artifacts'], receipt['source_hashes']]:
        for name, digest in mapping.items():
            bind(pins, name, digest)
    details = {}
    if args.stage == 'preparation':
        assert receipt['status'] == 'prepared_exact_qualified_full_scalar_v6_sampler_execution'
        assert (receipt['native_roles'], receipt['quartets'], receipt['effective_inputs'], receipt['aliases']) == (1620, 405, 135, 324)
        assert receipt['copied_job_file_matches_qualified_reference']
        assert receipt['expected_original_prerequisite_closures'] == 5
        assert not receipt['native_execution_launched'] and not receipt['posterior_qualified']
        assert sha(receipt['plan']) == receipt['plan_sha256']
        details.update(native_roles=1620, native_execution_launched=False, posterior_qualified=False)
    else:
        assert receipt['status'] == 'completed_protected_factor_probe_without_native_fault'
        assert receipt['debugger_exit_code'] == 0 and receipt['target_receipt_present']
        assert not receipt['native_sigsegv_observed'] and not receipt['fault_addresses']
        plan = json.loads(Path(receipt['plan']).read_text())
        assert sha(receipt['plan']) == receipt['plan_sha256']
        target_path = Path(plan['debugger_output']) / 'target_receipt.json'
        target = json.loads(target_path.read_text())
        assert target['status'] == 'completed_protected_factor_timing_probe_diagnostic_v7'
        assert target['planned_groups'] == target['executed_groups'] == 40
        assert target['probe_status_counts'] == {'timed_all_declared_points_agree_only': 40}
        assert target['inputs_preserved_for_all_attempted_groups'] and not target['detected_input_mutation']
        assert not target['unattempted_group_ids'] and not target['scientific_eligibility']
        for mapping in [target['source_hashes'], target['artifacts']]:
            for name, digest in mapping.items():
                bind(pins, name, digest)
        bind(pins, target_path)
        details.update(executed_groups=40, inputs_preserved=True,
            native_fault_observed=False, original_corruption_resolved=False,
            native_write_instruction_identified=False)
    inv = execution['invocation_id']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', args.unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [row for row in rows if inv in [row.get('_SYSTEMD_INVOCATION_ID'), row.get('USER_INVOCATION_ID')]]
    exact = [row for row in rows if row.get('_PID') == str(execution['wrapper']['pid'])
        and row.get('_CMDLINE') == ' '.join(execution['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(row['_SYSTEMD_INVOCATION_ID'] == inv for row in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=execution['wrapper'], invocation_id=inv)
    terminal = {key: value for key, value in execution.items()
        if key not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and 'Started ' in row.get('MESSAGE', '')]
    ends = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and row.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = args.execution.with_suffix('') / 'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + '\n')
    for path in [Path(__file__), args.execution, args.receipt, journal]:
        bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_' + args.stage.replace('-', '_') + '_execution_wait_zero',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_tool_session_id=args.original_tool_session, actual_tool_terminal_exit_code=0,
        invocation_id=inv, unit=args.unit, wrapper=execution['wrapper'],
        entire_wrapper_initial_and_terminal_payload_matched=True,
        original_start_records=1, original_completion_records=1,
        source_hashes=pins, **details, original_jobs_restarted=False, production_tolerance_changed=False,
        scientific_eligibility=False, biological_fits=0, gpu=False,
        scope='Original previously observed terminal-zero tool handle, full original wrapper initial/terminal '
              'journal payloads and manager start/end verified; tool transport ID is supplied from its '
              'earlier actual wait, not inferred from receipt. Source/artifact hashes freshly verified. '
              'Generic wrapper kernel label is superseded by its actual command. Protected run changes '
              'execution/allocation context and does not establish repair or the native write instruction. '
              'Post-exit GDB inspection errors are preserved. No original restart or biological acceptance.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__':
    main()
