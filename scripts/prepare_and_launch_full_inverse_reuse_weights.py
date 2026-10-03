#!/usr/bin/env python3
"""Gate and launch all original-cohort reuse controls, never weighted fits."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha
from full_inverse_reuse_weight_sources import SOURCES
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from run_full_inverse_reuse_weights import PRODUCER, READER, SUMMARY


def main():
    software = ['metadata/inverse_reuse_weight_software_' + k + '_20261003_v1.json'
                for k in ['validation','execution','transport']]
    v,e,t = [json.loads(Path(p).read_text()) for p in software]
    assert v['status'] == 'passed_complete_declared_inverse_reuse_controls_sql_fraction_contracts_v1'
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and not e['timed_out']
    assert e['receipt_sha256'] == sha(software[0])
    assert t['actual_tool_session_id'] == 90378 and t['actual_tool_terminal_exit_code'] == 0
    assert t['wrapper'] == e['wrapper'] and t['invocation_id'] == e['invocation_id']
    assert t['exact_wrapper_pid_journal_entries'] == 2
    assert t['original_start_records'] == t['original_completion_records'] == 1
    assert (v['mathematical_grouping_cases'],v['invalid_primitive_inputs_rejected'],
        v['synthetic_logical_cases'],v['synthetic_cohorts'],v['synthetic_case_row_occurrences']) == (40,3,24,5,40)
    assert len(v['rehashed_output_cases_rejected']) == 17 and len(v['rehashed_source_cases_rejected']) == 6
    assert v['completed_producer_reader_restarts_rejected'] is v['byte_exact_positive_restoration'] is True
    assert v['source_and_journal_fixtures_synthetic'] is True and v['scientific_eligibility'] is False
    for r in [v,e,t]: verify(r['source_hashes']);verify(r.get('artifacts',{}))
    rows = [json.loads(s) for s in subprocess.check_output(
        ['journalctl','--user','-u',t['unit'],'-o','json','--no-pager'],text=True).splitlines()]
    exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid'])
        and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])
        and r.get('_SYSTEMD_INVOCATION_ID') == e['invocation_id']]
    assert len(exact) == 2 and json.loads(exact[0]['MESSAGE'])['original_wrapper'] == e['wrapper']
    assert json.loads(exact[1]['MESSAGE'])['exit_code'] == 0
    parent_plans = dict(cases='metadata/full_matching_case_index_plan_20261002.json',
        covariance='metadata/full_expanded_covariance_plan_20261002.json',
        operator='metadata/full_entity_operator_plan_20261002.json',
        design='metadata/full_expanded_model_designs_plan_20261002_v2.json',
        cone='metadata/full_nonuniform_covariance_cone_plan_20261003_v1.json')
    parent_completions = dict(cases='metadata/full_matching_case_index_completed_20261002.json',
        covariance='metadata/full_expanded_covariance_completed_20261002.json',
        operator='metadata/full_entity_operator_bank_completed_20261002.json',
        design='metadata/full_expanded_model_designs_v2_completed_20261002.json',
        cone='metadata/full_nonuniform_covariance_cones_completed_20261003_v1.json')
    archives=[]
    for k,p in parent_completions.items():
        c=json.loads(Path(p).read_text())
        assert c['status']==SOURCES[k] and c['exact_process_journals_checked']==2
        assert c['scientific_eligibility'] is False and c['logical_cases']==75188
        assert sha(c['full_hash_archive'])==c['full_hash_archive_sha256'];archives.append(c['full_hash_archive'])
    design=json.loads(Path(parent_completions['design']).read_text())
    assert design['unique_cohorts']==4340 and design['cohort_member_occurrences']==34110120
    assert design['largest_cohort']==22934
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=16,
        swap_gib=0,blas_threads=1,address_space_gib=12,cpu_seconds_per_stage=21600,
        per_file_limit_mib=256,minimum_free_disk_gib=128,output_allowance_gib=8,
        original_logical_cases=75188,original_cohorts=4340,original_case_row_occurrences=34110120,
        control_policies=4,case_control_occurrences=136440480,
        raw_array_budget_bytes=34110120*96,largest_raw_cohort_array_bytes=22934*96,
        cohort_metadata_budget_bytes=4340*32768,
        runtime_planning_seconds_per_stage=[60,14400],runtime_calibrated=False,
        available_ram_gib=psutil.virtual_memory().available/2**30,
        free_disk_gib=shutil.disk_usage('.').free/2**30,new_cost_usd=0,gpu=False,
        scope='Full-data catalog generation and independent SQL/Fraction readback. Uncompressed array budget includes one row vector, three reuse-count vectors, four weights and four diagonals. Output allowance is a planning estimate, not a directory quota. CPU/AS/per-file/cgroup caps enforced; planning runtime is uncalibrated. No new structure inference or weighted fits.')
    assert psutil.virtual_memory().available>=16*2**30 and shutil.disk_usage('.').free>=128*2**30
    rp='metadata/full_inverse_reuse_weight_resources_20261003_v1.json';create(rp,resources)
    modules=['inverse_reuse_weight_controls','independent_inverse_reuse_weights','full_inverse_reuse_weight_sources',
        'run_full_inverse_reuse_weights','check_full_inverse_reuse_weights','run_inverse_reuse_weight_software_stage',
        'prepare_and_launch_full_inverse_reuse_weights','full_exact_covariance_sources',
        'full_expanded_model_design_sources','reference_measurement_union_sources','ancestral_chain_attempt',
        'launch_baliphy_reference_sampler_qualification','launch_full_triad_sequence_geometry_followup',
        'run_after_verified_dependencies_v2','close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    paths=['scripts/'+m+'.py' for m in modules]+[rp,*software,*parent_plans.values(),*parent_completions.values(),
        *archives,sys.executable,'/usr/bin/prlimit']
    output='results/phylogeny/full-original-cohort-inverse-reuse-controls-20261003-v1'
    assert not Path(output).exists()
    scope=('All original4340cohorts/75188logical cases/34110120cohort-row occurrences, uniform and three '
        'mean-one inverse-reuse sensitivity controls;136440480case/control occurrences. Cohort-local counts '
        'use original background nodes, versioned physical background pairs and connected family components. '
        'Preserve original QC membership, logical-case order and multiplicity links; do not expand selection '
        'records into observations. Reciprocal diagonal is a working residual-variance sensitivity assumption, '
        'not calibrated confidence-to-variance, literal power weighting, inverse-probability weighting, '
        'independent observation count or posterior ESS. Exact uniform-one classes retain a separate future '
        'qualification disposition. Fresh scoped parent archive/completion proofs and every consumed case '
        'table, caseID order, cohort membership and cone membership record are checked; unrelated archived '
        'X/y/operators/native artifacts are not freshly replayed. Every exported row/control checked using '
        'independent SQL partition counts and Fraction division. Two new actual original journals plus full '
        'stage source/artifact hashes gate closure. Actual weighted raw/REML numerical qualification, '
        'timing, fits and biological acceptance remain incomplete; all8aims remain open. No GPU/restart/new cost.')
    pp='metadata/full_inverse_reuse_weight_plan_20261003_v1.json'
    completion='metadata/full_inverse_reuse_weights_completed_20261003_v1.json'
    plan={k+'_plan':p for k,p in parent_plans.items()};plan.update({k+'_completion':p for k,p in parent_completions.items()})
    plan.update(expected=dict(logical_cases=75188,cohorts=4340,case_row_occurrences=34110120),output=output,
        resources=resources,pins={p:sha(p) for p in paths},completion=completion,scope=scope)
    create(pp,plan)
    prefix=['/usr/bin/prlimit','--as='+str(12*2**30),'--cpu=21600','--fsize='+str(256*2**20),'--']
    producer=launch('full-inverse-reuse-weights-v1',prefix+[sys.executable,'scripts/run_full_inverse_reuse_weights.py','--plan',pp],[],pp,cpus=2,memory=16)
    reader=launch('full-inverse-reuse-weights-v1-readback',prefix+[sys.executable,'scripts/run_full_inverse_reuse_weights.py','--plan',pp,'--reader'],[producer],pp,cpus=2,memory=16)
    cp='metadata/full_inverse_reuse_weight_completion_plan_20261003_v1.json'
    create(cp,dict(source_plan=pp,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER,reader_status=READER,completed_status='complete_verified_full_original_cohort_inverse_reuse_controls_v1',
        summary_fields=SUMMARY,launches=[producer,reader],
        pins={p:sha(p) for p in [pp,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=completion,scope=scope))
    closer=launch('full-inverse-reuse-weights-v1-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp,cpus=2,memory=16)
    create('metadata/full_inverse_reuse_weight_launches_20261003_v1.json',dict(status='queued_complete_original_cohort_reuse_control_export_and_independent_readback',
        source_plan=pp,source_plan_sha256=sha(pp),launches=[producer,reader,closer],
        expected=plan['expected'],expected_case_control_occurrences=136440480,
        residual_diagonal_prepared=False,raw_reml_basis_qualification_complete=False,
        weighted_fitting_launched=False,gpu=False,new_cost_usd=0,all_eight_aims_incomplete=True))


if __name__=='__main__':main()
