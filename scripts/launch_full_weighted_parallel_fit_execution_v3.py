#!/usr/bin/env python3
"""Prepare/submit complete original fit roles only after operational admission."""
import argparse
import json
from pathlib import Path
import sys
import shutil
import psutil

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from launch_baliphy_reference_sampler_qualification import launch
from prepare_full_weighted_parallel_fit_execution_v3 import run as prepare_operations, requirements
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify
from weighted_parallel_fit_execution_runtime_v3 import configuration

REQUEST='metadata/full_weighted_parallel_fit_execution_request_20261004_v3.json'
CONTROLLER='metadata/full_weighted_parallel_fit_controller_plan_20261004_v3.json'


def gate():
    gp=Path('metadata/full_weighted_parallel_fit_controller_software_validation_20261004_v3.json')
    ep=Path('metadata/full_weighted_parallel_fit_controller_software_execution_20261004_v3.json')
    tp=Path('metadata/full_weighted_parallel_fit_controller_software_transport_20261004_v3.json')
    v,e,t=[json.loads(p.read_text()) for p in [gp,ep,tp]]
    assert v['status']=='passed_full_weighted_parallel_fit_controller_scope_and_actual_native_boundary_contracts_v3'
    assert v['production_fitting_launched_or_queued'] is v['scientific_eligibility'] is False
    assert v['actual_native_boundary_probes']==4 and v['actual_successful_boundary_probes']==1 and v['actual_preserved_boundary_failures']==3
    assert v['completed_and_failed_attempt_restarts_refused'] is True and v['actual_missing_prerequisite_refusal'] is True
    assert v['successful_operational_preparation_branch_qualified'] is True
    assert v['successful_branch_closures_explicitly_mocked'] is True
    assert v['rejected_scope_and_capacity_changes']>=30 and v['rejected_native_custody_changes']>=12
    assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out'] and e['receipt_sha256']==sha(gp)
    assert t['actual_tool_terminal_exit_code']==0 and t['validation_sha256']==sha(gp)
    assert t['wrapper']==e['wrapper'] and t['invocation_id']==e['invocation_id']
    assert t['exact_wrapper_pid_journal_entries']==2 and t['original_start_records']==t['original_completion_records']==1
    for r in [v,e,t]:verify(r['source_hashes']);verify(r.get('artifacts',{}))
    return [str(gp),str(ep),str(tp)]


def prepare():
    request=json.loads(Path(REQUEST).read_text());verify(request['pins'])
    fit,timing,missing=requirements(request)
    if missing:return prepare_operations(REQUEST)
    proofs=gate()
    if not Path(request['operational_output']).exists():prepare_operations(REQUEST)
    operation=json.loads(Path(request['operational_output']).read_text());verify(operation['pins'])
    assert not Path(CONTROLLER).exists() and not Path(fit['output']).exists()
    modules=['weighted_parallel_fit_execution_runtime_v3','close_full_weighted_parallel_fit_execution_v3','launch_full_weighted_parallel_fit_execution_v3',
        'full_weighted_parallel_fit_admission_v1','prepare_full_weighted_parallel_fit_execution_v3','launch_baliphy_reference_sampler_qualification',
        'run_after_verified_dependencies_v2','close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    pins={p:sha(p) for p in [REQUEST,request['operational_output'],*operation['pins'],*proofs,
        *['scripts/'+m+'.py' for m in modules],sys.executable,'/usr/bin/prlimit']}
    plan=dict(status='software_qualified_full_weighted_fit_controller_not_launched',operation=request['operational_output'],
        expected=fit['expected'],original_closure_launches=request['original_closure_launches'],
        runtime_root='results/phylogeny/full-four-control-parallel-fit-execution-20261004-v3',
        completion='metadata/full_weighted_parallel_shared_entity_fits_completed_20261004_v3.json',
        launch_inventory='metadata/full_weighted_parallel_fit_controller_launches_20261004_v3.json',pins=pins,
        role_launch_inventory='metadata/full_weighted_parallel_fit_controller_role_launches_20261004_v3.json',
        controller_implemented=True,launch_state='not_launched_or_queued',scientific_eligibility=False,gpu=False,new_cost_usd=0,
        scope='All original20.832million candidates/49.7664million setting links with unchanged measured model plan. Full original numerical/timing/admission gates, exact role commands and native/cgroup/BLAS caps; every original nonfit/review/error retained. Separate producer/reader native identity/cap/log custody and all scientific source/output hashes plus original two controller journals close together. No automatic completed/failed attempt retry, skipped policy/method/tree or biological acceptance. Capacity is neither ETA nor convergence guarantee; inferential calibration/phylogenetic framework and all eight aims remain required.')
    assert not Path(plan['runtime_root']).exists();verify(pins);atomic(Path(CONTROLLER),plan)
    configuration(CONTROLLER)
    return plan


def queue():
    gate();plan,operation,fit,resources,bindings=configuration(CONTROLLER)
    for path in [fit['output'],plan['runtime_root'],plan['completion'],plan['launch_inventory'],plan['role_launch_inventory']]:assert not Path(path).exists()
    assert psutil.virtual_memory().available>=resources['minimum_available_ram_gib']*2**30
    assert shutil.disk_usage('.').free>=resources['minimum_free_disk_gib']*2**30
    for path in plan['original_closure_launches']:
        r=json.loads(Path(path).read_text());r['launch']=path;assert sha(r['plan'])==r['plan_sha256']
        assert fingerprint(r) is None; journal_terminal(r)
    producer=launch('full-weighted-parallel-fit-controller-v3',[sys.executable,'scripts/weighted_parallel_fit_execution_runtime_v3.py','--plan',CONTROLLER,'--role','producer'],
        plan['original_closure_launches'],CONTROLLER,cpus=2,memory=32)
    reader=launch('full-weighted-parallel-fit-controller-v3-readback',[sys.executable,'scripts/weighted_parallel_fit_execution_runtime_v3.py','--plan',CONTROLLER,'--role','reader'],
        [producer],CONTROLLER,cpus=2,memory=32)
    # This immutable two-role inventory exists before the closer can start.
    # The observer's three-controller inventory is written separately below.
    atomic(Path(plan['role_launch_inventory']),dict(controller_plan=CONTROLLER,controller_plan_sha256=sha(CONTROLLER),launches=[producer,reader]))
    closer=launch('full-weighted-parallel-fit-controller-v3-closure',[sys.executable,'scripts/close_full_weighted_parallel_fit_execution_v3.py','--plan',CONTROLLER],
        [producer,reader],CONTROLLER,cpus=2,memory=32)
    atomic(Path(plan['launch_inventory']),dict(status='queued_full_weighted_fits_after_verified_original_closures',
        controller_plan=CONTROLLER,controller_plan_sha256=sha(CONTROLLER),launches=[producer,reader,closer],
        fit_plan=operation['fit_plan'],fit_plan_sha256=operation['fit_plan_sha256'],expected=plan['expected'],
        actual_fit_resources_installed=True,scientific_eligibility=False,gpu=False,new_cost_usd=0))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--action',choices=['prepare','queue'],required=True);a=p.parse_args()
    result=prepare() if a.action=='prepare' else queue()
    if result is not None:print(json.dumps(result,indent=2),flush=True)
