#!/usr/bin/env python3
"""Recheck captured factor changes and corrected scalar fixtures without a restart."""
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/current_analysis_error_recheck_20261004_user_0748.json')
    assert not output.exists()
    diagnostic_path = Path('metadata/retained_timing_ordered_probe_diagnostic_20261004_v2.json')
    transport_path = Path('metadata/retained_timing_ordered_probe_diagnostic_transport_20261004_v2.json')
    diagnostic = json.loads(diagnostic_path.read_text())
    transport = json.loads(transport_path.read_text())
    assert transport['validation_sha256'] == sha(diagnostic_path)
    assert transport['actual_tool_terminal_exit_code'] == 0
    assert transport['entire_terminal_payload_matched']
    pins = dict(transport['source_hashes'])
    verify(pins)
    captures = diagnostic['detected_input_mutation']['factor_value_change_capture']
    assert set(captures) == {'pmsf_profile_profile'}
    capture = captures['pmsf_profile_profile']
    observed_path = Path(capture['observed_snapshot'])
    baseline_path = Path('results/phylogeny/retained-timing-covariance-diagnostic-20261004-v1/pmsf_profile_profile-factor.npy')
    assert sha(observed_path) == capture['observed_snapshot_sha256']
    baseline = np.load(baseline_path, mmap_mode='r', allow_pickle=False)
    observed = np.load(observed_path, mmap_mode='r', allow_pickle=False)
    assert baseline.dtype == observed.dtype == np.dtype('float64')
    assert baseline.shape == observed.shape == (22881, 301)
    changed = np.argwhere(baseline.view(np.uint64) != observed.view(np.uint64))
    assert len(changed) == capture['changed_entries'] == 8
    assert changed.tolist() == [[8513, column] for column in range(115, 123)]
    for entry, coordinate in zip(capture['samples'], changed):
        row, column = map(int, coordinate)
        assert entry['row'] == row and entry['column'] == column
        assert baseline[row, column].item().hex() == entry['baseline_hex']
        assert observed[row, column].item().hex() == entry['observed_hex']
    assert not capture['current_probe_uses_changed_tree']
    qualified_path = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    qualified = json.loads(qualified_path.read_text())
    verify(qualified['source_hashes'])
    comparisons = [compare_tsv(item['new_directory']) for item in qualified['paired_prior_checks']]
    assert sum(item['rows'] for item in comparisons) == 63
    assert sum(item['mapped_values_compared'] for item in comparisons) == 2709
    failure_path = Path('metadata/full_retained_shared_entity_timing_failure_20261004_v2.json')
    failure = json.loads(failure_path.read_text())
    verify(failure['source_hashes'])
    root = Path('results/phylogeny/full-retained-shared-entity-input-timing-20261003-v1')
    assert len(list((root / 'cohorts').glob('*.receipt.json'))) == 10
    assert not (root / 'receipt.json').exists()
    startup_path = Path('metadata/baliphy_scalar_v6_preflight_execution_checkpoint_20261004_user_0744.json')
    startup = json.loads(startup_path.read_text())
    for path in [Path(__file__), diagnostic_path, transport_path, baseline_path,
                 observed_path, qualified_path, failure_path, startup_path]:
        bind(pins, path)
    verify(pins)
    result = dict(status='completed_user_saved_error_evidence_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        retained_timing_original_failure_unresolved=True,
        original_retained_timing_partial_cohorts=10,
        captured_changed_factor_entries_rechecked=8,
        changed_tree='pmsf_profile_profile', probe_tree=capture['current_probe_tree'],
        changed_coordinates=changed.tolist(),
        affected_snapshot=str(observed_path), affected_snapshot_sha256=sha(observed_path),
        diagnostic_terminal_exit_code=0, diagnostic_executed_groups=diagnostic['executed_groups'],
        diagnostic_unattempted_groups=len(diagnostic['unattempted_group_ids']),
        corrected_v6_saved_rows_rechecked=63, corrected_v6_mapped_values_rechecked=2709,
        corrected_v6_comparisons=comparisons,
        startup_checked_utc=startup['checked_utc'],
        startup_unclosed_passing_roles=startup['unclosed_startup_checkpoints'],
        startup_expected_roles=1620,
        source_hashes=pins, root_cause_established=False, original_jobs_restarted=False,
        new_native_runs=0, new_mcmc_runs=0, gpu=False, biological_fits=0,
        scope='Fresh read-only source/artifact hash checks and independent bitwise comparison '
              'of saved diagnostic factor against unchanged baseline. Captured changes revalidated; '
              'original timing failure remains unresolved. Corrected retained scalar fixtures '
              'reread independently. This is not another native replay, a root-cause diagnosis, '
              'a repair, full source-archive replay, full startup closure or biological acceptance.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ['source_hashes', 'corrected_v6_comparisons']}, indent=2))


if __name__ == '__main__':
    main()
