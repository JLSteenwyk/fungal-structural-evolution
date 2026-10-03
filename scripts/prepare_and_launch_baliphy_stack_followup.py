#!/usr/bin/env python3
"""Prepare scoped correction after full source closures and software execution."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha
from baliphy_stack_followup import build_jobs, prerequisites, SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from run_baliphy_stack_followup import PRODUCER_STATUS, READER_STATUS, COMPLETED_STATUS


def main():
    gate_path = Path('metadata/baliphy_stack_followup_software_validation_20261003_v1.json')
    gate = json.loads(gate_path.read_text()); verify(gate['source_hashes'])
    assert gate['status'] == 'passed_full_role_stack_followup_software_contracts'
    assert (gate['full_original_roles'],gate['followup_roles'],gate['fresh_disjoint_seeds']) == (1620,24,24)
    assert len(gate['design_cases_rejected']) == 9 and len(gate['serialization_and_resource_cases_rejected']) == 10
    assert gate['full_producer_reader_serialization_passed'] is gate['native_execution_mocked'] is True
    # Terminal tool session3980 exited zero. Bind the actual original native
    # software controller journal; its synthetic fixtures remain explicit.
    unit = 'fungal-baliphy-stack-followup-software-20261003-v1.service'
    invocation = '8adf514103724f3db55187d689149f2d'
    raw = subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
    jr = [json.loads(line) for line in raw.splitlines()]
    software_root = Path('data/software_audits/baliphy-stack-followup-software-20261003-v1')
    controller = json.loads((software_root/'pipeline/controller.json').read_text())
    exact = [r for r in jr if r.get('_PID') == str(controller['pid']) and r.get('_CMDLINE') == ' '.join(controller['cmdline'])]
    assert exact and {r['_SYSTEMD_INVOCATION_ID'] for r in exact} == {invocation}
    resources = [r for r in jr if r.get('USER_INVOCATION_ID') == invocation and r.get('CPU_USAGE_NSEC')]
    assert resources and not any('Failed with result' in r.get('MESSAGE','') or 'Main process exited' in r.get('MESSAGE','') for r in jr if r.get('USER_INVOCATION_ID') == invocation)
    journal = software_root/'original-invocation-journal.jsonl'; assert not journal.exists(); journal.write_text(raw)
    execution_path = 'metadata/baliphy_stack_followup_software_execution_20261003_v1.json'
    create(execution_path,dict(status='verified_original_stack_followup_software_tool_wait_exited_zero',
        session_id=3980,actual_tool_terminal_exit_code=0,invocation_id=invocation,
        original_controller={k:controller[k] for k in ['pid','created','cmdline']},original_process_messages=len(exact),
        original_completion_resource_records=len(resources),software_receipt=str(gate_path),software_receipt_sha256=sha(gate_path),
        source_hashes={str(gate_path):sha(gate_path),str(journal):sha(journal)},
        configured_resources=dict(cpus=2,memory_gib=8,swap_gib=0,address_space_gib=6,cpu_seconds=480,
            service_wall_seconds=600,per_file_limit_gib=1,blas_threads=1),
        scope='Actual original software tool wait/journal success. The pipeline cgroup/telemetry/native data were explicitly mocked; the recorded software launch caps are configuration, not retrospective live procfs measurements. Manager completion peak is not treated as native peak.'))
    source_path = 'metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json'
    observer_path = 'metadata/baliphy_joint_sampler_resource_observation_v3_plan_20261003.json'
    source = json.loads(Path(source_path).read_text())
    originals = json.loads(Path(source['jobs']).read_text())
    rows_path = Path(source['output'])/'dispositions.json'; rows = json.loads(rows_path.read_text())
    hp = Path(source['historical_sampler_plan']); historical_plan = json.loads(hp.read_text())
    historical = json.loads(Path(historical_plan['jobs']).read_text())
    diagnosis_path = 'metadata/baliphy_stack_followup_diagnosis_20261003_v1.json'
    diagnosis = json.loads(Path(diagnosis_path).read_text()); verify(diagnosis['source_hashes'])
    jobs = build_jobs(originals,rows,historical,diagnosis)
    inputs = Path('results/ancestral/baliphy-stack-followup-inputs-20261003-v1');inputs.mkdir(exist_ok=False)
    job_path = inputs/'jobs.json';create(str(job_path),jobs)
    for j in jobs:verify(j['config']['pins'])
    dependencies = [json.loads(Path(p).read_text())['launches'][2] for p in [source['launch_inventory'],
        'metadata/baliphy_joint_sampler_resource_observation_v3_launches_20261003.json']]
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=4,memory_gib=200,workers=4,
        reservation_capacity_gib=192,swap_gib=0,blas_threads=1,iterations=20,native_stack_bytes=64*2**20,
        native_address_space_gib=48,native_file_limit_gib=2,reader_cpus=2,reader_memory_gib=32,
        sum_native_wall_timeout_worker_hours=sum(j['config']['timeout_seconds'] for j in jobs)/3600,
        maximum_parallel_native_timeouts_upper_wall_hours=sum(j['config']['timeout_seconds'] for j in jobs)/3600/4,
        output_allowance_gib=32,minimum_free_disk_gib=256,maximum_followup_joint_frames=72,
        fresh_followup_roles=24,full_original_roles=1620,model_prior_input_unchanged=True,
        runtime_uncalibrated=True,finish_eta=None,planning_output_is_not_hard_quota=True,
        new_cost_usd=0,gpu=False,posterior_qualified=False,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        caveat='Initial-state corrections do not calibrate full20iteration runtime. Timeout sum is an upper allowance, not an ETA; dividing by four is an idealized balanced scheduling allowance, not a strict wall bound.192GiB leases remain through output audit under200GiB group cap. Only native stack setting and fresh seed differ.32GiB output planning excludes the original referenced dataset and is not a disk quota. No longer posterior or new infrastructure.')
    own = ['baliphy_stack_followup','run_baliphy_stack_followup','check_baliphy_stack_followup',
        'prepare_and_launch_baliphy_stack_followup','validate_baliphy_stack_followup_diagnosis',
        'ancestral_chain_attempt','baliphy_joint_sampler_qualification_v3','baliphy_joint_sampler_gates_v3',
        'baliphy_sampler_resource_observation','run_baliphy_joint_sampler_qualification_v3','reference_sampler_memory_budget',
        'independent_joint_ancestral_frames','independent_native_ancestral_alignment','independent_native_ancestral_topology',
        'independent_short_sampler_outputs_v2','independent_short_sampler_outputs_v3','baliphy_reference_sampler_qualification',
        'readback_independent_baliphy_chain','run_after_verified_dependencies_v2','launch_baliphy_reference_sampler_qualification',
        'close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    paths = [Path('scripts/'+n+'.py') for n in own] + [gate_path,job_path,rows_path,Path(source_path),Path(observer_path),
        Path(diagnosis_path),Path(execution_path),Path(source['jobs']),hp,Path(historical_plan['jobs']),Path(source['mapping']),*map(Path,dependencies)]
    pins = {str(p):sha(p) for p in paths}
    scope = ('Complete1620original roles/405quartets/135inputs/324aliases retained with explicit original and follow-up provenance. '
        'Twenty-four new20iteration attempts cover all two-input/three-prior/four-chain diagnosedSIGSEGV roles;64MiB stack is '
        'scoped to their native prlimit commands. Fresh seed namespace disjoint all original/current/historical/known diagnostic '
        'seeds; no continuation,concatenation,restart,replacement of originals or automatic repeated attempts. Original models, '
        'priors,alignments,trees,AS48GiB,file2GiB and CPU/wall caps preserved. Complete source sampler/resource/artifact/two-journal '
        'closures and two causal strict initial-frame diagnostics required. All1596original successes reused by immutable '
        'reference; all24original failures retained even if follow-ups succeed. All follow-up native/config/identity/scalar/'
        'joint-frame/array checks, sampled stack/memory/cgroup limits, FIFO reservations and independent serialized full-role '
        'readback required. Failed fresh attempts and unresolved quartets remain explicit. No adequate posterior, precise final '
        'native resource peaks, long-run qualification, historical bad_alloc repair, accepted tree/root/model or completed '
        'biological aims implied. GPU prediction remains paused; no paid resources or global defaults changed.')
    plan_path = 'metadata/baliphy_stack_followup_plan_20261003_v1.json'
    output = 'results/ancestral/full-baliphy-stack-followup-20261003-v1'
    plan = dict(source_plan=source_path,source_resource_plan=observer_path,diagnosis=diagnosis_path,
        jobs=str(job_path),mapping=source['mapping'],output=output,pins=pins,resources=resources,dependencies=dependencies,
        software_validation=str(gate_path),completion='metadata/baliphy_stack_followup_completed_20261003_v1.json',
        launch_inventory='metadata/baliphy_stack_followup_launches_20261003_v1.json',scope=scope)
    # Full source archives, original sampler failures and resource observations
    # are checked before writing the immutable plan or starting new inference.
    prerequisites(plan)
    assert psutil.virtual_memory().available >= 200*2**30 and shutil.disk_usage('.').free >= 256*2**30
    assert not Path(output).exists();create(plan_path,plan)
    producer = launch('baliphy-stack-followup-v1',[sys.executable,'scripts/run_baliphy_stack_followup.py','--plan',plan_path],
        dependencies,plan_path,cpus=4,memory=200)
    reader = launch('baliphy-stack-followup-v1-readback',[sys.executable,'scripts/run_baliphy_stack_followup.py','--plan',plan_path,'--reader'],
        [producer],plan_path,cpus=2,memory=32)
    cp = 'metadata/baliphy_stack_followup_completion_plan_20261003_v1.json'
    create(cp,dict(source_plan=plan_path,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,completed_status=COMPLETED_STATUS,
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={p:sha(p) for p in [plan_path,producer,reader,'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=scope))
    closer = launch('baliphy-stack-followup-v1-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],
        [producer,reader],cp,cpus=2,memory=32)
    create(plan['launch_inventory'],dict(status='queued_full_role_scoped_stack_followup',source_plan=plan_path,
        source_plan_sha256=sha(plan_path),launches=[producer,reader,closer],original_roles=1620,fresh_followup_roles=24,
        historical_failures_replaced=False,longer_posterior_launched=False,posterior_qualified=False,gpu=False,new_cost_usd=0))
    print(json.dumps(resources),flush=True)


if __name__ == '__main__':main()
