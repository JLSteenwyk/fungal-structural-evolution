#!/usr/bin/env python3
"""Full-scope prospective fitting arithmetic; no launch or accepted input claim."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    paths=['metadata/full_expanded_model_designs_plan_20261002_v2.json',
        'metadata/full_uniform_covariance_qualification_plan_20261002.json',
        'metadata/full_expanded_model_inputs_completed_20261002.json',
        'metadata/full_entity_operator_bank_completed_20261002.json',
        'metadata/full_uniform_covariance_qualification_resources_20261002.json',
        'metadata/shared_entity_optimizer_validation_20261002_v2.json',
        'metadata/independent_shared_entity_likelihood_validation_20261002_v5.json',
        'scripts/shared_entity_likelihood.py','scripts/fit_shared_entity_likelihood.py',
        'scripts/independent_shared_entity_likelihood.py','environments/shared-entity-likelihood-20261002.yml',
        'metadata/full_expanded_covariance_completed_20261002.json',
        'scripts/prepare_full_shared_entity_fit_resources.py']
    bindings={s:sha(s) for s in paths}
    plan=json.loads(Path(paths[0]).read_text());qualification=json.loads(Path(paths[1]).read_text())
    receipt_path=Path(plan['output'])/'receipt.json';bindings[str(receipt_path)]=sha(receipt_path)
    design=json.loads(receipt_path.read_text())
    assert design['plan_sha256']==sha(paths[0]) and design['scientific_eligibility'] is False
    assert design['status']=='complete_full_expanded_model_designs_pending_independent_readback'
    assert design['model_setting_rows']==qualification['expected']['model_setting_rows']==622080
    assert design['unique_cohorts']==4340 and design['unique_designs']==130200 and design['unique_fit_inputs']==260400
    for path,expected in [(paths[2],'complete_verified_full_expanded_model_inputs'),(paths[3],'complete_verified_full_entity_operator_bank')]:
        closed=json.loads(Path(path).read_text());assert closed['status']==expected
        assert sha(closed['full_hash_archive'])==closed['full_hash_archive_sha256']
        bindings[closed['full_hash_archive']]=closed['full_hash_archive_sha256']
    for path in paths[5:7]:
        validation=json.loads(Path(path).read_text());assert validation['status'].startswith('passed_')
        for source,digest in validation['source_hashes'].items():assert sha(source)==digest;bindings[source]=digest
    modes=2;trees=len(design['trees']);methods=['ml','reml'];outcomes=2
    assert trees==5 and design['unique_fit_inputs']==design['unique_designs']*outcomes
    unique_candidates=design['unique_fit_inputs']*modes*trees*len(methods)
    setting_links=design['model_setting_rows']*modes*trees*len(methods)
    starts=3;max_evaluations=500;ratios=6
    numerical=json.loads(Path(paths[4]).read_text())
    full_species_one=numerical['maximum_five_full_species_factor_bytes']//5
    # Peak allocations need actual full-stage measurement. These array counts
    # are conservative ceilings within one cohort, not a process RSS bound.
    covariance=json.loads(Path('metadata/full_expanded_covariance_completed_20261002.json').read_text())
    assert covariance['status']=='complete_verified_full_expanded_covariance' and covariance['trees']==design['trees']
    column_batch=32;n=design['largest_cohort'];factor_rank=covariance['rank'];assert factor_rank==301
    bytes_per_candidate=64*1024;bytes_per_audit=32*1024;bytes_per_link=1024
    record_bound=unique_candidates*(bytes_per_candidate+bytes_per_audit)+setting_links*bytes_per_link
    result=dict(status='prospective_full_shared_entity_fit_resource_inventory_pending_closed_qualification_and_timing',
        checked_utc=datetime.now(timezone.utc).isoformat(),full_original_fixed_settings=design['model_setting_rows'],
        observed_original_cohorts=design['unique_cohorts'],observed_original_designs=design['unique_designs'],
        observed_original_fit_inputs=design['unique_fit_inputs'],observed_scope_not_yet_independently_closed=True,
        methods=methods,loading_modes=['signed','unsigned'],working_trees=design['trees'],outcomes=outcomes,
        maximum_unique_uniform_candidates=unique_candidates,complete_setting_method_mode_tree_links=setting_links,
        possible_covariance_variance_ratios=ratios,deterministic_optimizer_starts=starts,
        per_start_search_evaluation_budget=max_evaluations,per_start_final_replay_budget=1,
        maximum_search_and_final_replay_evaluations=unique_candidates*starts*(max_evaluations+1),
        maximum_candidate_replay_and_local_curvature_evaluations=unique_candidates*(2+4*ratios),
        additional_independent_global_optimization_budget_not_yet_specified=True,
        eligible_fit_count_pending_full_uniform_covariance_qualification=True,
        proposed_cpu_quota_per_worker=2,proposed_memory_gib_per_worker=32,proposed_swap_gib=0,proposed_blas_threads=1,
        largest_original_cohort=n,species_factor_rank=factor_rank,
        maximum_single_component_matrix_bytes=numerical['maximum_separate_component_bytes'],
        maximum_five_cached_entity_kernel_bytes=numerical['maximum_five_kernel_component_stack_bytes'],
        maximum_full_bank_one_species_factor_bytes=full_species_one,
        maximum_cohort_species_factor_bytes=n*factor_rank*8,
        independent_entity_column_batch=column_batch,maximum_dense_streamed_kernel_image_bytes=n*column_batch*8,
        candidate_record_budget_bytes=bytes_per_candidate,independent_audit_record_budget_bytes=bytes_per_audit,
        setting_link_budget_bytes=bytes_per_link,conservative_uncompressed_record_budget_bytes=record_bound,
        conservative_uncompressed_record_budget_gib=record_bound/2**30,proposed_output_scratch_reserve_gib=1024,
        free_disk_gib_at_observation=shutil.disk_usage(Path.cwd()).free/2**30,
        runtime_calibrated=False,runtime_eta_hours=None,production_launch_authorized_by_this_inventory=False,
        gpu_usage=False,new_paid_resources=False,source_hashes=bindings,
        scope='Complete prospective uniform ML/REML fitting grid, including all original setting links and '
            'exact shared-input reuse; qualification/review rows must remain even when no optimizer runs. '
            'Counts are preliminary producer scope, not accepted designs or fitted effects. Evaluation '
            'budgets are maxima, not expected iterations or wall time. Output estimates require record '
            'size monitoring; array ceilings are not process peak-memory guarantees. Before production '
            'launch close original design/qualification stages, finish restartable producer/independent '
            'reader/export contracts, estimate full-scope timing with qualified inputs, specify independent '
            'global-optimization checks and enforce resource/disk limits. Nonuniform/control variants, '
            'uncertainty calibration, accepted phylogenetic framework and all eight aims remain required.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['status','maximum_unique_uniform_candidates',
        'complete_setting_method_mode_tree_links','maximum_search_and_final_replay_evaluations',
        'conservative_uncompressed_record_budget_gib','runtime_calibrated']}))


if __name__=='__main__':main()
