#!/usr/bin/env python3
"""Inventory and queue the complete timing census behind original numeric closure."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import psutil

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from full_weighted_covariance_sources_v2 import POLICIES, MODES
from full_expanded_model_design_sources import DEGREES, AXES
from full_expanded_model_input_sources import ORDERS, NUISANCE
from launch_baliphy_reference_sampler_qualification import launch
from prepare_and_launch_full_weighted_covariance_qualification import software
from prepare_full_weighted_shared_entity_timing_v2 import PRODUCER, READER
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify

FIT='metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json'
PLAN='metadata/full_weighted_shared_entity_timing_plan_20261004_v2.json'
RESOURCES='metadata/full_weighted_shared_entity_timing_resources_20261004_v2.json'
DEPENDENCY='metadata/full_weighted_covariance_numerical_v1_closure_launch_20261003.json'
SUMMARY=['logical_cases','model_setting_rows','unique_cohorts','candidate_rows','eligible_candidates','timing_groups',
    'source_status_counts','timing_status_counts','measured_candidate_coverage','unmeasured_review_candidate_coverage',
    'conditional_budget_weighted_seconds']


def gates():
    value,paths=software('full_weighted_timing','20261004_v3',50869)
    assert value['status']=='passed_complete_four_control_timing_census_native_numeric_and_accounting_contracts_v3'
    assert (value['total_candidate_rows'],value['full_grid_selected_groups'],value['actual_additional_numeric_groups'],
        value['actual_additional_numeric_points'])==(72000,320,64,128)
    assert value['full_grid_native_probes_mocked'] is value['scientific_eligibility'] is False
    assert value['source_and_journal_fixtures_synthetic'] is True and value['fits_computed']==0
    assert value['primary_evaluation_budget']==1503 and value['independent_evaluation_budgets']==[1522,1526,1530]
    assert len(value['rehashed_output_alterations_rejected'])==15 and len(value['malformed_probe_exports_rejected'])==15
    assert len(value['synthetic_source_closure_alterations_rejected'])==12
    for k in ['closed_producer_chunks_replayed_without_rewrite','completed_producer_and_reader_restarts_refused',
        'full_grid_selections_independently_recomputed','review_input_arrays_content_addressed_and_byte_replayed',
        'review_input_array_corruption_rejected']:assert value[k] is True
    return paths


def prepare():
    proofs=gates();fp=Path(FIT);fit=json.loads(fp.read_text());verify(fit['pins'])
    assert fit['launch_state']=='not_launched_or_queued' and not Path(fit['output']).exists()
    assert fit['optimizer']['gradient_tolerance']==1e-6 and fit['independent_audit']['column_batch']==32
    assert fit['policies']==POLICIES and fit['loading_modes']==MODES
    qp=Path(fit['qualification_plan']);q=json.loads(qp.read_text());verify(q['pins'])
    expected=q['expected'];assert expected['case_row_occurrences']==34110120
    dp=Path(q['design_plan']);design=json.loads((Path(json.loads(dp.read_text())['output'])/'receipt.json').read_text())
    cp=Path('metadata/full_expanded_covariance_completed_20261002.json');cov=json.loads(cp.read_text())
    assert cov['rank']==301 and cov['trees']==fit['trees']
    assert design['unique_cohorts']==4340 and design['unique_designs']==130200 and design['largest_cohort']==22934
    assert design['cohort_member_occurrences']==34110120 and design['trees']==fit['trees']
    designs_per_cohort=len(ORDERS)*len(AXES)*len(DEGREES);columns=1+max(DEGREES)+len(NUISANCE)
    assert (designs_per_cohort,columns)==(30,10)
    groups=expected['cohorts']*4*2*5*2*2;assert groups==694400
    candidates=expected['designs']*2*4*2*5*2;assert candidates==fit['expected']['candidate_rows']==20832000
    # Content-addressed original arrays avoid multiplying the source matrices
    # by methods, modes and responses in every review. No compression ratio or
    # cross-cohort sharing is assumed for these planning byte allowances.
    numeric_columns=designs_per_cohort*columns+designs_per_cohort*2+len(POLICIES)+len(fit['trees'])*cov['rank']
    arrays=expected['case_row_occurrences']*(numeric_columns*8+8+64)
    record_bytes=candidates*2048+groups*32768*2
    budget=arrays+record_bytes
    reserve=768;minimum=868;assert budget<reserve*2**30
    resources=dict(status='prospective_complete_weighted_timing_resources',prepared_utc=datetime.now(timezone.utc).isoformat(),
        cpus=2,memory_gib=32,swap_gib=0,blas_threads=1,address_space_gib=24,native_cpu_seconds_per_stage=604800,
        per_file_limit_gib=2,minimum_free_disk_gib=minimum,output_scratch_allowance_gib=reserve,
        maximum_original_candidates=candidates,maximum_eligible_groups=groups,scaled_variance_points=[0.,1.],
        maximum_primary_evaluations_producer=groups*2,maximum_independent_evaluations_producer=groups*2,
        maximum_primary_evaluations_reader=groups*2,maximum_independent_evaluations_reader=groups*2,
        all_groups_fresh_numeric_replay_in_reader=True,review_inputs_content_addressed=True,
        original_case_row_occurrences=expected['case_row_occurrences'],designs_per_cohort=designs_per_cohort,
        maximum_active_columns=columns,species_factor_columns=cov['rank'],maximum_original_cohort=design['largest_cohort'],
        maximum_review_array_raw_bytes=arrays,census_and_probe_and_review_record_budget_bytes=record_bytes,
        total_uncompressed_planning_gib=budget/2**30,maximum_original_single_factor_bytes=design['largest_cohort']*cov['rank']*8,
        primary_fit_evaluation_budget=1503,independent_fit_evaluation_budgets={'4':1522,'5':1526,'6':1530},
        runtime_planning_core_hours_per_stage=[6,168],runtime_calibrated=False,production_finish_eta=None,
        native_cpu_cap_is_eta=False,wall_time_limit_enforced=False,available_ram_gib=psutil.virtual_memory().available/2**30,
        free_disk_gib=shutil.disk_usage('.').free/2**30,gpu=False,new_cost_usd=0,fits_computed=0,
        scope='Whole20.832millioncandidate census, maximum694400eligible groups and both numeric points, actualfour control diagonals; every group independently recomputed during numeric/accounting readback. Content-addressed review arrays byte-bound, no compression or cross-cohort sharing assumed. Worst all-review matrix allowance uses all30designs/60responses/4diagonals/5rank301factors plus rows/labels at each original cohort.768GiB scratch is planning allowance, not enforced directory quota; free-disk checks and native CPU/AS/cgroup/per-file/BLAS limits apply.168CPUhour cap is a capacity allocation, not measured timing or finish ETA. Source loading/export/hash/review-array I/O not inferred from point measurements. No native fit, new charge, GPU or biological acceptance.')
    assert psutil.virtual_memory().available>=64*2**30 and shutil.disk_usage('.').free>=minimum*2**30
    atomic(Path(RESOURCES),resources)
    modules=['full_weighted_shared_entity_timing','full_weighted_timing_contracts_v2','prepare_full_weighted_shared_entity_timing_v2',
        'readback_full_weighted_shared_entity_timing_v2','check_full_weighted_shared_entity_timing_v3','check_full_weighted_shared_entity_timing_v2',
        'prepare_and_launch_full_weighted_shared_entity_timing','record_weighted_timing_software_transport',
        'launch_baliphy_reference_sampler_qualification','run_after_verified_dependencies_v2','close_full_triad_sequence_stage',
        'record_completed_process_handoffs_v2','full_weighted_shared_entity_fit_sources','weighted_shared_entity_candidate',
        'independent_positive_diagonal_basis_context','shared_entity_likelihood','independent_shared_entity_likelihood',
        'independent_shared_entity_likelihood_fast','readback_full_covariance_qualification','full_weighted_fit_exports']
    paths=['scripts/'+m+'.py' for m in modules]+proofs+[FIT,RESOURCES,str(qp),str(dp),str(cp),DEPENDENCY,sys.executable,'/usr/bin/prlimit']
    scope='Complete4340cohorts/130200designs/260400response inputs and20.832millioncandidate identities, allfourpolicies/twomodes/fivetrees/twooutcomes/MLandREML. Original exclusions/reviews retained; deterministic maximum active-column/condition/candidate-ID representative per eligiblecohort/control/mode/tree/outcome/method. Freshactual-Dcomponent versus independentlatent qualification and primary versus independentcomponent-spectral probes atscaledvariance0and1. All groups receive fresh numerical replay in reader; hardware durations are observed, not independently reproducible. Exact original numeric closer2975938 gates execution; all consumed sources/output bytes and original two journals gate timing closure. Conditional budgets include source validation, both production guards and reader qualification, optimizer/failure-replay paths, but neither bound unknown optimizer-point costs nor forecast finish time. All construction/precision reviews retain original arrays through shared content-addressed custody. No missing-control/method/tree reduction, pilot, native optimization, fit resource installation, prior-job restart, GPU, paid resources or accepted biological effect. Full calibration/phylogenetic framework and all eight aims remain required.'
    plan=dict(status='software_qualified_complete_weighted_input_timing_pending_original_numerical_closure',
        fit_plan=FIT,fit_plan_sha256=sha(FIT),scaled_variance_points=[0.,1.],numerical_closure_launch=DEPENDENCY,
        output='results/phylogeny/full-four-control-shared-entity-input-timing-20261004-v2',
        completion='metadata/full_weighted_shared_entity_timing_completed_20261004_v2.json',
        closure_plan='metadata/full_weighted_shared_entity_timing_completion_plan_20261004_v2.json',
        launch_inventory='metadata/full_weighted_shared_entity_timing_launches_20261004_v2.json',
        resources=resources,summary_fields=SUMMARY,pins={p:sha(p) for p in paths},scope=scope,
        prepared_utc=datetime.now(timezone.utc).isoformat(),production_fitting_launched=False,gpu=False,new_cost_usd=0)
    assert not Path(plan['output']).exists();verify(plan['pins']);atomic(Path(PLAN),plan)
    print(json.dumps(dict(status=plan['status'],candidate_rows=candidates,maximum_groups=groups,
        total_uncompressed_planning_gib=resources['total_uncompressed_planning_gib'],timing_plan=PLAN),indent=2))


def queue():
    gates();plan=json.loads(Path(PLAN).read_text());verify(plan['pins']);fit=json.loads(Path(FIT).read_text());verify(fit['pins'])
    assert sha(FIT)==plan['fit_plan_sha256'] and fit['launch_state']=='not_launched_or_queued' and not Path(fit['output']).exists()
    for path in [plan['output'],plan['completion'],plan['closure_plan'],plan['launch_inventory']]:assert not Path(path).exists()
    assert psutil.virtual_memory().available>=64*2**30 and shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    handle=json.loads(Path(DEPENDENCY).read_text());handle['launch']=DEPENDENCY
    assert sha(handle['plan'])==handle['plan_sha256']
    if fingerprint(handle) is None:journal_terminal(handle)
    prefix=['/usr/bin/prlimit','--as='+str(24*2**30),'--cpu=604800','--fsize='+str(2*2**30),'--']
    producer=launch('full-weighted-shared-entity-timing-v2',prefix+[sys.executable,
        'scripts/prepare_full_weighted_shared_entity_timing_v2.py','--plan',PLAN],[DEPENDENCY],PLAN,cpus=2,memory=32)
    root=Path(plan['output'])
    reader=launch('full-weighted-shared-entity-timing-v2-readback',prefix+[sys.executable,
        'scripts/readback_full_weighted_shared_entity_timing_v2.py','--plan',PLAN,'--output',str(root/'readback.json')],
        [producer],PLAN,cpus=2,memory=32)
    atomic(Path(plan['closure_plan']),dict(source_plan=PLAN,producer_receipt=str(root/'receipt.json'),
        independent_readback=str(root/'readback.json'),producer_status=PRODUCER,reader_status=READER,
        completed_status='complete_verified_full_four_control_input_timing_accounting_v2',summary_fields=SUMMARY,
        launches=[producer,reader],pins={p:sha(p) for p in [PLAN,producer,reader,'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=plan['scope']))
    closer=launch('full-weighted-shared-entity-timing-v2-closure',[sys.executable,
        'scripts/close_full_triad_sequence_stage.py','--plan',plan['closure_plan']],[producer,reader],plan['closure_plan'],cpus=2,memory=32)
    atomic(Path(plan['launch_inventory']),dict(status='queued_complete_weighted_input_timing_after_original_numerical_closure',
        source_plan=PLAN,source_plan_sha256=sha(PLAN),launches=[producer,reader,closer],candidate_rows=20832000,
        maximum_eligible_timing_groups=694400,production_fitting_launched=False,gpu=False,new_cost_usd=0,all_eight_aims_incomplete=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--action',choices=['prepare','queue'],required=True)
    prepare() if p.parse_args().action=='prepare' else queue()
