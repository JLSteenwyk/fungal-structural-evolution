#!/usr/bin/env python3
"""Record fresh formatter checks and the original completed V9 diagnostic."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind, verify


def close_execution(path, unit, tool):
    execution = json.loads(path.read_text())
    assert execution['exit_code'] == 0 and not execution['timed_out']
    assert tool['terminal']['exit_code'] == 0
    assert sha(execution['receipt']) == execution['receipt_sha256']
    verify(execution['source_hashes'])
    verify(execution['artifacts'])
    inv = execution['invocation_id']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [r for r in rows if inv in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID'))]
    wrapper = execution['wrapper']
    exact = [r for r in rows if r.get('_PID') == str(wrapper['pid'])
             and r.get('_CMDLINE') == ' '.join(wrapper['cmdline'])]
    assert len(exact) == 2
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=wrapper, invocation_id=inv)
    terminal = {k: v for k, v in execution.items()
                if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    return execution, dict(invocation_id=inv, wrapper=wrapper,
        original_tool_session_id=tool['original_tool_session_id'],
        original_tool_terminal_exit_code=0, whole_wrapper_payloads_matched=True,
        manager_start_records=1, manager_completion_records=1)


def main():
    output = Path('metadata/current_analysis_error_recheck_20261004_user_1122.json')
    assert not output.exists()
    payload_path = Path('metadata/current_analysis_error_original_tool_payloads_20261004_user_1122.json')
    payload = json.loads(payload_path.read_text())
    fresh_path = Path('metadata/baliphy_logging_error_recheck_execution_20261004_user_1122.json')
    diagnostic_path = Path('metadata/retained_factor_watched_traced_debugger_execution_20261004_v1.json')
    fresh, fresh_proof = close_execution(fresh_path,
        'fungal-baliphy-logging-error-recheck-20261004-user-1122.service', payload['fresh_check'])
    diagnostic, diagnostic_proof = close_execution(diagnostic_path,
        'fungal-retained-factor-watched-traced-debugger-20261004-v1.service', payload['diagnostic'])
    probe = json.loads(Path(fresh['receipt']).read_text())
    original = probe['numeric_probes']['original']
    assert original['error_reproduced'] and len(original['incorrect_constant_indices']) == 5
    assert probe['numeric_probes']['cjson']['native_exact_roundtrips'] == 12
    qualified_path = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    qualified = json.loads(qualified_path.read_text())
    verify(qualified['source_hashes'])
    comparisons = [compare_tsv(item['new_directory']) for item in qualified['paired_prior_checks']]
    assert sum(c['rows'] for c in comparisons) == 63
    assert sum(c['mapped_values_compared'] for c in comparisons) == 2709
    target_path = Path('results/phylogeny/retained-factor-watched-traced-debugger-20261004-v1/target_receipt.json')
    target = json.loads(target_path.read_text())
    debugger = json.loads(Path(diagnostic['receipt']).read_text())
    assert target['executed_groups'] == target['planned_groups'] == 40
    assert not target['unattempted_group_ids'] and target['detected_input_mutation'] is None
    assert target['inputs_preserved_for_all_attempted_groups']
    assert target['probe_status_counts'] == {'timed_all_declared_points_agree_only': 40}
    assert debugger['debugger_exit_code'] == 0 and not debugger['native_write_hit']
    verify(target['source_hashes'])
    failure_path = Path('metadata/full_retained_shared_entity_timing_failure_20261004_v2.json')
    failure = json.loads(failure_path.read_text())
    verify(failure['source_hashes'])
    failed_root = Path('results/phylogeny/full-retained-shared-entity-input-timing-20261003-v1')
    assert len(list((failed_root / 'cohorts').glob('*.receipt.json'))) == 10
    assert not (failed_root / 'receipt.json').exists()
    pins = {}
    for path in [Path(__file__), payload_path, fresh_path, diagnostic_path, Path(fresh['receipt']),
                 Path(diagnostic['receipt']), qualified_path, target_path, failure_path]:
        bind(pins, path)
    result = dict(status='completed_fresh_user_error_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        original_formatter_error_reproduced=True, incorrect_native_constants=5,
        example=dict(expected=2.34e-10, observed=original['values'][0]),
        corrected_native_roundtrips=12, corrected_v6_saved_rows_rechecked=63,
        corrected_v6_mapped_values_rechecked=2709, comparisons=comparisons,
        diagnostic_completed_groups=40, diagnostic_total_groups=40,
        diagnostic_native_write_hit=False, diagnostic_inputs_preserved=True,
        retained_timing_original_failure_unresolved=True, original_failed_partial_cohorts=10,
        memory_corruption_resolved=False, root_cause_established=False,
        fresh_probe_execution_proof=fresh_proof, diagnostic_execution_proof=diagnostic_proof,
        source_hashes=pins, new_mcmc_runs=0, original_jobs_restarted=False,
        installed_software_changed=False, production_tolerance_changed=False,
        biological_fits=0, scientific_eligibility=False, gpu=False,
        scope='The installed formatter still fails the fresh native constant test. The corrected '
              'project path passes twelve constants and independent readback of 2709 retained values. '
              'The original V9 diagnostic completed all forty groups without detected corruption; '
              'instrumentation changes execution context and this is not evidence of a root-cause repair. '
              'The original failed timing run remains preserved and incomplete.')
    verify(pins)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'comparisons']}, indent=2))


if __name__ == '__main__':
    main()
