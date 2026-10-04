#!/usr/bin/env python3
"""Fresh full archive rehash and original terminal checks for sampler/telemetry."""
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/baliphy_scalar_v6_sampler_full_closure_recheck_20261004_v1.json')
    assert not output.exists()
    completions = []
    handles = []
    pins = {}
    for label, path in [('sampler', Path('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json')),
                        ('observer', Path('metadata/baliphy_scalar_v6_resource_observation_completed_20261004_v1.json'))]:
        value = json.loads(path.read_text())
        archive_path = Path(value['full_hash_archive'])
        assert sha(archive_path) == value['full_hash_archive_sha256']
        archive = json.loads(archive_path.read_text())
        assert len(archive['source_hashes']) == value['bound_source_hashes']
        assert len(archive['services']) == value['exact_process_journals_checked'] == 2
        assert value['scientific_eligibility'] is archive['scientific_eligibility'] is False
        verify(archive['source_hashes'])
        for service in archive['services']:
            observed = observe(service['launch'])
            assert observed['status'] == 'verified_original_terminal_success_with_bound_completed_artifacts'
            handles.append(observed)
        for key in ['producer_receipt', 'independent_readback']:
            assert sha(value[key]) == value[key+'_sha256']
        for p in [path, archive_path]:bind(pins, p)
        completions.append(dict(stage=label, completion=str(path), sha256=sha(path),
            archive=str(archive_path), archive_sha256=sha(archive_path), bindings_rehashed=value['bound_source_hashes'],
            original_journals=2, status=value['status']))
    main_completion = json.loads(Path(completions[0]['completion']).read_text())
    assert main_completion['full_chains'] == 1620 and main_completion['full_quartets'] == 405
    assert main_completion['status_counts'] == {
        'full_joint_short_sampler_output_integrity_checked_not_posterior':1584,
        'explicit_nonfinite_or_literal_null_scalar_output_retained_for_review':12,
        'unsuccessful_sampler_qualification_attempt_retained':24}
    root = Path('results/ancestral/full-baliphy-scalar-v6-short-sampler-20261004-v1')
    failures = []
    for path in sorted((root/'chains').glob('*.json')):
        row = json.loads(path.read_text())
        if row['status'] == 'unsuccessful_sampler_qualification_attempt_retained':
            assert row['exit_code'] == -11 and not row['scalar_integrity_accepted']
            failures.append(dict(chain_id=row['chain_id'], family=row['family'], prior=row['prior_label'],
                role=row['chain_role'], effective_input_group=row['effective_input_group'], exit_code=row['exit_code'],
                native_receipt=row['native_receipt'], native_receipt_sha256=row['native_receipt_sha256']))
            bind(pins, path)
    assert len(failures) == 24
    bind(pins, Path(__file__))
    result = dict(status='fresh_full_scalar_v6_sampler_and_resource_archives_verified',
        checked_utc=datetime.now(timezone.utc).isoformat(), closures=completions,
        original_terminal_handles=handles, full_roles=1620, full_quartets=405,
        finite_integrity_roles=1584, special_value_reviews=12, native_failed_roles=24,
        failure_family_counts=dict(Counter(r['family'] for r in failures)), failed_roles=failures,
        scalar_mapped_values_checked=main_completion['scalar_v6_mapped_values_checked'],
        joint_saved_frames=main_completion['joint_saved_frames'],
        joint_ancestral_residue_category_pairs=main_completion['joint_ancestral_residue_category_pairs'],
        source_hashes=pins, posterior_qualified=False, scientific_eligibility=False, gpu=False,
        original_attempts_restarted=False,
        scope='All116217sampler and6501telemetry archive bindings freshly rehashed; original '
              'four producer/reader terminal journals rechecked by exact original handles. '
              'All24failed roles remain retained, separate from twelve special-value reviews. '
              'The full original independent reader already replayed native output integrity; '
              'this is a fresh whole-archive/terminal check, not another numeric/native sampler '
              'execution, a crash repair or adequate posterior qualification.')
    verify(pins)
    with output.open('x') as handle:json.dump(result, handle, indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_terminal_handles','failed_roles']}, indent=2))


if __name__ == '__main__':main()
