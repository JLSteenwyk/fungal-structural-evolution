#!/usr/bin/env python3
"""Retain fresh native probes and distinguish correction from unresolved failure."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import psutil

from ancestral_chain_attempt import sha
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/current_analysis_error_recheck_20261004_user_0650.json')
    assert not output.exists()
    execution = Path('metadata/baliphy_logging_error_recheck_execution_20261004_user_0650.json')
    e = json.loads(execution.read_text())
    assert e['exit_code'] == 0 and not e['timed_out']
    assert e['status'] == 'exited_zero_with_receipt'
    receipt = Path(e['receipt'])
    assert sha(receipt) == e['receipt_sha256']
    probe = json.loads(receipt.read_text())
    assert probe['numeric_probes']['original']['error_reproduced']
    assert len(probe['numeric_probes']['original']['incorrect_constant_indices']) == 5
    assert probe['numeric_probes']['cjson']['native_exact_roundtrips'] == 12
    pins = {}
    for mapping in [e['source_hashes'], e['artifacts'], probe['source_hashes']]:
        for path, digest in mapping.items():
            bind(pins, path, digest)
    unit = 'fungal-logging-error-recheck-20261004-user-0650.service'
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    inv = e['invocation_id']
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID')]]
    exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid'])
             and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(r['_SYSTEMD_INVOCATION_ID'] == inv for r in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: v for k, v in e.items()
                if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = execution.with_suffix('') / 'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + '\n')
    qualified = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    corrected = json.loads(qualified.read_text())
    assert corrected['status'] == 'passed_scalar_json_v6_full_source_and_paired_native_qualification'
    for path, digest in corrected['source_hashes'].items():
        bind(pins, path, digest)
    comparisons = [compare_tsv(item['new_directory']) for item in corrected['paired_prior_checks']]
    assert sum(c['rows'] for c in comparisons) == 63
    assert sum(c['mapped_values_compared'] for c in comparisons) == 2709
    adapter_path = Path('metadata/baliphy_scalar_v6_joint_adapter_software_transport_20261004_v2.json')
    adapter = json.loads(adapter_path.read_text())
    assert adapter['status'] == 'verified_original_v6_scalar_joint_adapter_v2_software_wait_zero'
    for path, digest in adapter['source_hashes'].items():
        bind(pins, path, digest)
    failed = Path('metadata/full_retained_shared_entity_timing_failure_20261004_v2.json')
    failure = json.loads(failed.read_text())
    for path, digest in failure['source_hashes'].items():
        bind(pins, path, digest)
    launches = [failure['original_producer_launch']] + [
        c['original_launch'] for c in failure['current_dependency_controllers']]
    for launch in launches:
        try:
            process = psutil.Process(launch['pid'])
            assert process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            pass
    root = Path('results/phylogeny/full-retained-shared-entity-input-timing-20261003-v1')
    assert len(list((root / 'cohorts').glob('*.receipt.json'))) == 10
    assert not (root / 'receipt.json').exists()
    diagnostic_path = Path('metadata/retained_timing_covariance_diagnostic_checkpoint_20261004_user_0649.json')
    diagnostic = json.loads(diagnostic_path.read_text())
    weighted_path = Path('metadata/full_weighted_covariance_qualification_execution_checkpoint_20261004_user_0650.json')
    weighted = json.loads(weighted_path.read_text())
    for path in [Path(__file__), Path('scripts/read_baliphy_scalar_json_v6b.py'),
                 execution, receipt, journal, qualified, adapter_path, failed,
                 diagnostic_path, weighted_path]:
        bind(pins, path)
    verify(pins)
    result = dict(status='completed_fresh_user_analysis_error_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        original_formatter_error_reproduced=True, incorrect_native_constants=5,
        original_small_value_example=dict(expected=2.34e-10,
            observed=probe['numeric_probes']['original']['values'][0]),
        cjson_twelve_native_roundtrips=12, corrected_v6_saved_rows_rechecked=63,
        corrected_v6_mapped_values_rechecked=2709, corrected_v6_comparisons=comparisons,
        v6_admission_adapter_source_and_artifact_bindings_rechecked=len(adapter['source_hashes']),
        full_grid_v6_native_execution_launched=False,
        retained_timing_original_failure_unresolved=True, retained_timing_partial_cohorts=10,
        retained_timing_failure_reexecuted=False,
        diagnostic_checked_utc=diagnostic['checked_utc'],
        diagnostic_last_progress=diagnostic['last_progress_lines'],
        diagnostic_terminal=False, diagnostic_receipt_present=diagnostic['receipt_present'],
        weighted_checked_utc=weighted['checked_utc'],
        weighted_last_progress=weighted['last_progress_message'],
        weighted_failure_files=weighted['current_failure_files'],
        actual_tool_session_id=10189, actual_tool_terminal_exit_code=0,
        invocation_id=inv, wrapper=e['wrapper'], entire_terminal_payload_matched=True,
        original_start_records=1, original_completion_records=1,
        source_hashes=pins, gpu=False, new_mcmc_runs=0, original_jobs_restarted=False,
        scope='Fresh original native formatter/CJSON probes and independent retained V6 scalar '
              'readback. V6 qualification and adapter source/artifact pins freshly rehashed; '
              'full-grid V6 startup/controller/closure gates remain incomplete. Original retained '
              'timing failure evidence remains unchanged with no completion; it was not rerun. '
              'Diagnostic/weighted snapshots refer to their explicitly dated live observers. '
              'Generic wrapper kernel scope is superseded by its actual logging-recheck command. '
              'No installed software repair, covariance acceptance or biological inference.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__':
    main()
