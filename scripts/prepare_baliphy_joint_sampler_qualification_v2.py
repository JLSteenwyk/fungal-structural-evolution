#!/usr/bin/env python3
"""Prepare all 1620 joint-logger short qualification roles; no native launch."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil

import psutil

from ancestral_chain_attempt import sha,write_json
from baliphy_joint_sampler_qualification_v2 import build_jobs
from run_baliphy_reference_preflight import verify


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True); p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args(); gate=json.loads(a.validation.read_text())
    assert gate['status']=='passed_full_joint_sampler_qualification_v2_software_contracts'
    assert gate['full_roles']==1620 and gate['full_quartets']==405 and gate['full_source_configurations_checked']
    assert gate['actual_retained_native_fixture_roles']==3 and gate['actual_retained_native_fixture_frames']==9
    assert gate['full_mock_producer_reader_serialization_checked'] and gate['artificial_failed_roles_retained']==2
    assert {'missing_frame_export','missing_export_directory'} <= set(gate['serialization_cases_rejected'])
    assert len(gate['altered_designs_rejected'])==11 and len(gate['prerequisite_gate_negatives'])==7
    verify({'pins':gate['source_hashes']})
    future_path=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future=json.loads(future_path.read_text());verify(future)
    native_path=Path(future['source_plan']);native=json.loads(native_path.read_text());verify(native)
    originals=json.loads(Path(native['jobs']).read_text()); roles=json.loads(Path(future['future_roles']).read_text())
    jobs=build_jobs(roles,originals); root=a.inputs.resolve();root.mkdir(exist_ok=False)
    job_path=root/'jobs.json';write_json(job_path,jobs)
    dependencies=['metadata/baliphy_joint_logger_preflight_closure_launch_20261003.json',
        'metadata/baliphy_reference_sampler_qualification_closure_launch_20261003.json',
        'metadata/baliphy_sampler_resource_observation_closure_launch_20261003.json',
        'metadata/independent_short_sampler_replay_v2_closure_launch_20261003.json']
    stage_plans=['metadata/baliphy_joint_logger_preflight_plan_20261003.json',str(native_path),
        'metadata/baliphy_sampler_resource_observation_plan_20261003.json',
        'metadata/independent_short_sampler_replay_v2_plan_20261003.json']
    owned=['baliphy_joint_sampler_qualification_v2','baliphy_joint_sampler_gates',
        'run_baliphy_joint_sampler_qualification_v2','check_baliphy_joint_sampler_qualification_v2',
        'prepare_baliphy_joint_sampler_qualification_v2','launch_baliphy_joint_sampler_qualification_v2',
        'record_baliphy_joint_sampler_qualification_v2_checkpoint','independent_joint_ancestral_frames',
        'reference_sampler_memory_budget','independent_native_ancestral_topology',
        'independent_native_ancestral_alignment','independent_short_sampler_outputs_v2',
        'baliphy_reference_sampler_qualification','readback_independent_baliphy_chain','ancestral_chain_attempt']
    paths=[Path('scripts',name+'.py') for name in owned]+[a.validation,future_path,Path(future['future_roles']),
        native_path,Path(native['jobs']),Path(native['mapping']),job_path,*map(Path,dependencies),*map(Path,stage_plans)]
    pins=dict(future['pins']);pins.update({str(path):sha(path) for path in paths});verify({'pins':pins})
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),workers=16,cpus=16,memory_gib=200,
        reservation_capacity_gib=192,swap_gib=0,blas_threads=1,
        role_address_space_gib_counts={'12':1512,'48':108},risk_families=native['resources']['risk_families'],
        iterations=20,per_file_limit_gib=2,output_allowance_gib=128,minimum_free_disk_gib=256,
        maximum_joint_frames=4860,maximum_candidate_frames=19440,
        candidate_state_category_array_component_bytes=587929824,
        historical_old_initialization_linear_twenty_iteration_worker_seconds=native['resources']['linear_old_successful_worker_seconds'],
        sum_wall_timeout_worker_hours=sum(j['config']['timeout_seconds'] for j in jobs)/3600,
        runtime_uncalibrated=True,finish_eta=None,new_cost_usd=0,gpu=False,posterior_qualified=False,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage(root).free/2**30,
        caveat='Full20iteration computational gate, not adequate posterior. Old1000iteration linear scale has different initialization and is not an ETA; configured timeouts are upper allowances. AS leases include native execution/output audit but are declared reservations, not measured Python/native peaks. EightGiB group headroom, enforced200GiB/no swap and disk guard;128GiB output is planning, not a hard global quota. Candidate arrays exclude raw all-node JSON, coordinates, free residues and legacy logs. Waiting for four closed prerequisites avoids simultaneous old/new native grids.')
    scope=('All1620preparedfuturejointroles/405quartets/135inputs/324aliases:20iterations, three joint frames per successful role. '
        'Startup seeds used for firstMCMC, not continued posterior samples. All four original prerequisite closures mandatory, '
        'including corrected ambiguity reader and observed resource boundaries. Same-record all-node states/categories, rooted '
        'clades/candidates, observed-tip constraints, anchored/unanchored residue coordinates, complete native/config/process/artifact '
        'and exact arrays checked. Failed/invalid rows and unresolved quartets retained without automatic retry. Readback never '
        'recreates missing exports; producer binding verification required. Complete source/artifact/two-original-journal closure '
        'mandatory. No longer posterior, mixing/likelihood/model/root/predictor acceptance, historical allocation repair, GPU, '
        'charges or completed biological aims implied. Existing sources/attempts/caps unchanged.')
    plan=dict(jobs=str(job_path),output='results/ancestral/full-baliphy-joint-sampler-qualification-20261003-v2',
        mapping=native['mapping'],startup_plan=stage_plans[0],historical_sampler_plan=stage_plans[1],
        historical_resource_plan=stage_plans[2],historical_replay_plan=stage_plans[3],dependencies=dependencies,
        completion='metadata/baliphy_joint_sampler_qualification_v2_completed_20261003.json',
        completion_plan='metadata/baliphy_joint_sampler_qualification_v2_completion_plan_20261003.json',
        launch_inventory='metadata/baliphy_joint_sampler_qualification_v2_launches_20261003.json',
        software_validation=str(a.validation),pins=pins,resources=resources,scope=scope)
    with a.plan.open('x') as handle:handle.write(json.dumps(plan,indent=2)+'\n')
    print(json.dumps(dict(plan=str(a.plan),**resources)),flush=True)


if __name__=='__main__':main()
