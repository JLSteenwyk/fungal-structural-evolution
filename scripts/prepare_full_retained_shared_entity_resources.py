#!/usr/bin/env python3
"""Inventory the complete retained fitting grid before timing or fitting launch."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify


def inventory():
    paths = [
        'metadata/full_expanded_model_designs_plan_20261002_v2.json',
        'metadata/full_uniform_covariance_qualification_resources_20261002.json',
        'metadata/full_expanded_covariance_completed_20261002.json',
        'metadata/full_reduced_covariance_qualification_launches_20261003_v3.json',
        'metadata/full_retained_fitting_software_validation_20261003_v2.json',
        'metadata/full_retained_fitting_software_transport_20261003_v2.json',
        'metadata/full_retained_timing_software_validation_20261003_v1.json',
        'metadata/full_retained_timing_software_transport_20261003_v1.json',
        'scripts/prepare_full_retained_shared_entity_resources.py',
    ]
    plan = json.loads(Path(paths[0]).read_text())
    receipt_path = Path(plan['output']) / 'receipt.json'
    design = json.loads(receipt_path.read_text())
    assert design['plan_sha256'] == sha(paths[0])
    assert (design['unique_cohorts'], design['unique_designs'], design['unique_fit_inputs'],
            design['model_setting_rows']) == (4340, 130200, 260400, 622080)
    arrays = json.loads(Path(paths[1]).read_text())
    covariance = json.loads(Path(paths[2]).read_text())
    assert covariance['rank'] == 301 and covariance['trees'] == design['trees']
    retained = json.loads(Path(paths[3]).read_text())
    assert retained['real_certificate_basis_counts'] == {'4': 7232, '5': 1448}
    bindings = {p: sha(p) for p in paths + [str(receipt_path)]}
    for p in paths[4:8]:
        gate = json.loads(Path(p).read_text())
        verify(gate['source_hashes'])
        bindings.update(gate['source_hashes'])
    candidates = design['unique_fit_inputs'] * 2 * 5 * 2
    links = design['model_setting_rows'] * 2 * 5 * 2
    groups = design['unique_cohorts'] * 2 * 5 * 2 * 2
    assert (candidates, links, groups) == (5208000, 12441600, 173600)
    record_bytes = candidates * (64 * 1024 + 32 * 1024) + links * 1024
    return dict(
        status='prospective_complete_retained_fitting_and_timing_resource_inventory',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        original_cohorts=4340, original_designs=130200, original_fit_inputs=260400,
        original_fixed_settings=622080, potential_candidates=candidates, setting_fit_links=links,
        methods=['ml', 'reml'], loading_modes=['signed', 'unsigned'], trees=design['trees'],
        exact_certificate_basis_counts=retained['real_certificate_basis_counts'],
        retained_covariance_components=[4, 5], variance_ratio_dimensions=[3, 4],
        primary_evaluations_per_candidate=1503,
        maximum_primary_evaluations=candidates * 1503,
        independent_evaluations_per_candidate_by_basis={'4': 1522, '5': 1526},
        maximum_independent_evaluations=candidates * 1526,
        candidate_record_budget_bytes=64 * 1024, independent_record_budget_bytes=32 * 1024,
        setting_link_budget_bytes=1024, uncompressed_record_budget_bytes=record_bytes,
        uncompressed_record_budget_gib=record_bytes / 2**30,
        fitting_output_scratch_reserve_gib=1024, fitting_minimum_free_disk_gib=1124,
        largest_original_cohort=design['largest_cohort'], species_factor_rank=301,
        maximum_component_matrix_bytes=arrays['maximum_separate_component_bytes'],
        maximum_cached_component_stack_bytes=arrays['maximum_five_kernel_component_stack_bytes'],
        maximum_cohort_species_factor_bytes=design['largest_cohort'] * 301 * 8,
        proposed_fitting_cpus_per_worker=2, proposed_fitting_memory_gib_per_worker=32,
        timing=dict(cpus_per_stage=2, memory_gib=32, swap_gib=0, blas_threads=1,
                    address_space_gib=24, native_cpu_seconds_per_stage=604800,
                    per_file_limit_gib=2, minimum_free_disk_gib=164,
                    output_allowance_gib=64, maximum_groups=groups,
                    scaled_variance_points=[0.0, 1.0],
                    maximum_initial_likelihood_evaluations=groups * 4,
                    maximum_primary_constructor_calls=groups,
                    maximum_independent_constructor_calls=groups,
                    maximum_full_source_validation_calls=groups,
                    maximum_fresh_backend_guard_calls=groups * 2,
                    maximum_latent_reader_guard_calls=groups,
                    claimed_review_groups_reproduced_by_reader=True),
        free_disk_gib=shutil.disk_usage('.').free / 2**30,
        available_ram_gib=psutil.virtual_memory().available / 2**30,
        source_loading_and_timing_not_yet_completed=True,
        runtime_calibrated=False, expected_runtime_hours=None, finish_eta=None,
        numerical_threads=1, new_cost_usd=0, gpu=False, production_fitting_launched=False,
        source_hashes=bindings,
        scope='Complete original uniform grid remains explicit, including all exclusions and reviews. '
              'Exact certificate counts are closed operator proofs, not accepted designs. Retained fitting '
              'keeps three starts and500search evaluations plus final replay, independent curvature and '
              'three spectral searches. Record allowances and component array ceilings are planning '
              'counts, not process RSS or enforced total output quotas. Timing uses two CPU/32GiB/no swap '
              'with24GiB native address space and168CPU-hour per-stage allocation; this limit is a '
              'capacity budget, not a runtime estimate or guaranteed finish. Constructor/source/both '
              'production guards and independent latent qualification are separately priced. Full source '
              'closure and whole-grid timing remain required before fitting. No global optimality proof, '
              'biological acceptance, pilot reduction, new paid resources or GPU prediction.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inventory()
    with args.output.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'potential_candidates', 'setting_fit_links',
                                         'uncompressed_record_budget_gib', 'runtime_calibrated']}))


if __name__ == '__main__':
    main()
