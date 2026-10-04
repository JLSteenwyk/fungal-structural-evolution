#!/usr/bin/env python3
"""Fresh full source/archive verification of completed paired predictor draws."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from matched_predictor_resampling import SUMMARY_FIELDS
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/matched_predictor_resampling_full_closure_verified_20261004_v1.json')
    assert not output.exists()
    cp = Path('metadata/matched_predictor_resampling_completed_20261003_v1.json')
    completed = json.loads(cp.read_text()); ap = Path(completed['full_hash_archive'])
    assert completed['status'] == 'complete_verified_full_matched_predictor_paired_resampling'
    assert sha(ap) == completed['full_hash_archive_sha256']
    archive = json.loads(ap.read_text())
    assert len(archive['source_hashes']) == completed['bound_source_hashes'] == 176656
    assert len(archive['services']) == completed['exact_process_journals_checked'] == 2
    verify(archive['source_hashes'])
    summary = {k: completed[k] for k in SUMMARY_FIELDS}
    assert summary == archive['summary']
    assert summary['unique_inputs'] == 133 and summary['resampling_cases'] == 53200
    assert summary['native_roles'] == 372400 and summary['unresolved_cases'] == 0
    assert summary['complete_cases'] == 53200
    assert summary['finite_branch_values'] == summary['serialized_branch_value_slots'] == 6146000
    assert summary['native_status_counts'] == {'native_point_output_integrity_checked_not_model_qualified': 372400}
    root = Path(completed['producer_receipt']).parent
    producer, reader = [json.loads(Path(completed[k]).read_text()) for k in ['producer_receipt', 'independent_readback']]
    for k in ['producer_receipt', 'independent_readback']:
        assert sha(completed[k]) == completed[k + '_sha256']
    assert reader['producer_receipt_sha256'] == sha(completed['producer_receipt'])
    assert all(producer[k] == reader[k] == v for k, v in summary.items())
    assert producer['scientific_eligibility'] is reader['scientific_eligibility'] is False
    manifest = json.loads((root / 'case_manifest.json').read_text())
    keys = sorted({r['original_input_id'] for r in manifest}); assert len(keys) == 133
    assert [(r['original_input_id'], r['mode'], r['replicate']) for r in manifest] == [
        (key, mode, i) for key in keys for mode in sorted(summary['modes']) for i in range(200)]
    assert len({r['case_id'] for r in manifest}) == 53200
    pins = {}; handles = []
    for index, service in enumerate(archive['services']):
        handle = observe(service['launch'])
        assert handle['status'] == 'verified_original_terminal_success_with_bound_completed_artifacts'
        launch = json.loads(Path(service['launch']).read_text())
        rows = [json.loads(line) for line in subprocess.check_output(
            ['journalctl', '--user', '-u', launch['unit'], '--all', '-o', 'json', '--no-pager'], text=True).splitlines()]
        matches = []
        for row in rows:
            try: value = json.loads(row['MESSAGE'])
            except (KeyError, ValueError, TypeError): continue
            if value == summary: matches.append(row)
        assert len(matches) == 1
        wrapper = [r for r in rows if r.get('_PID') == str(launch['pid'])
                   and r.get('_CMDLINE') == ' '.join(launch['cmdline'])]
        assert wrapper and {r['_SYSTEMD_INVOCATION_ID'] for r in wrapper} == {matches[0]['_SYSTEMD_INVOCATION_ID']}
        journal = Path(f'metadata/matched_predictor_resampling_{index}_original_whole_journal_20261004_v1.jsonl')
        with journal.open('x') as f:
            for row in rows: f.write(json.dumps(row, sort_keys=True) + '\n')
        handles.append(dict(original_handle=handle, native_summary_pid=matches[0]['_PID'],
                            native_summary_cmdline=matches[0]['_CMDLINE'], whole_native_summary_matched=True))
        for path in [Path(service['launch']), journal]: bind(pins, path)
    for path in [Path(__file__), cp, ap]: bind(pins, path)
    verify(pins)
    result = dict(status='fresh_full_matched_predictor_resampling_archive_and_original_journals_verified',
                  checked_utc=datetime.now(timezone.utc).isoformat(), **summary,
                  archive_bindings_freshly_rehashed=176656, original_terminal_handles=handles,
                  source_hashes=pins, paired_uncertainty_stage_numerically_closed=True,
                  predictor_overlap_markers=71, predictor_overlap_fungal_taxa=21,
                  representative_of_full_fungal_sampling=False,
                  calibrated_confidence_intervals=False, scientific_eligibility=False,
                  all_eight_aims_incomplete=True, original_jobs_restarted=False, gpu=False,
                  scope='All176656 declared source/output/archive bindings rehashed and original '
                        'wrapper/manager success evidence plus complete native terminal summaries '
                        'matched. Original reader independently reconstructed all53200 draws and '
                        '372400 archived native roles. Fresh verification is not another native fit '
                        'or archive-member numeric replay. Paired distribution summaries, Monte Carlo '
                        'precision, coverage/model assessment, direct structural controls and broad '
                        'fungal replication remain required before accepted effects or intervals.')
    with output.open('x') as f: json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_terminal_handles']}, indent=2))


if __name__ == '__main__':
    main()
