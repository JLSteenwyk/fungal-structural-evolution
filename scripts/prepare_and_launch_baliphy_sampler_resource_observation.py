#!/usr/bin/env python3
"""Start a separate read-only observer of the original full sampler queue."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from baliphy_sampler_resource_observation import SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint
from run_baliphy_reference_preflight import verify
from run_baliphy_sampler_resource_observation import PRODUCER_STATUS, READER_STATUS, COMPLETED_STATUS


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--validation',type=Path,required=True)
    args=parser.parse_args();gate=json.loads(args.validation.read_text());verify(dict(pins=gate['source_hashes']))
    assert gate['status']=='passed_full_sampler_resource_observation_software_contracts'
    assert gate['full_roles_checked']==1620 and gate['actual_live_native_observations']>0
    assert gate['full_producer_and_reader_serialization_checked']
    assert gate['full_missing_observation_accounting_checked'] and gate['stable_unrelated_command_rejected']
    assert gate['empty_or_transitioning_argv_not_accepted_as_native'] and gate['wrapper_not_mistaken_for_native']
    assert len(gate['malformed_or_altered_observations_rejected'])==9 and len(gate['serialization_rejection_cases'])==8
    native_path=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native=json.loads(native_path.read_text());verify(native)
    original_path='metadata/baliphy_reference_sampler_qualification_launch_20261003.json'
    original=json.loads(Path(original_path).read_text());assert fingerprint(original) is not None
    assert psutil.virtual_memory().available>=2*2**30 and shutil.disk_usage('.').free>=64*2**30
    scope=('Read-only full1620role resource observation for the exact original sampler controller and native attempts. '
        'PID/create/command/configuration/cgroup/limits bound; empty/transitioning argv, wrapper phases, process races '
        'and missed fast attempts remain unavailable observations, not accepted native readings or terminal proofs. '
        'Every role and failed/unstarted attempt retained. Linux-reported virtual/resident/CPU and cgroup memory/pressure '
        'snapshots are approximate and can miss later native peaks; no precise final per-attempt peak or long-chain '
        'memory guarantee. Observer never writes GPU/native cgroup limits or restarts/signals inference. '
        'Full streaming serialized replay/artifact/source/two original observer completion journals required for '
        'observational closure, not sampler/model/tree/root/posterior or biological acceptance. All eight aims incomplete.')
    own=['baliphy_sampler_resource_observation','run_baliphy_sampler_resource_observation',
        'check_baliphy_sampler_resource_observation','check_baliphy_sampler_resource_serialization',
        'prepare_and_launch_baliphy_sampler_resource_observation','record_baliphy_sampler_resource_observation_checkpoint']
    paths=[Path('scripts/'+name+'.py') for name in own]+[args.validation,native_path,Path(original_path)]
    pins=dict(native['pins']);pins.update({str(p):sha(p) for p in paths})
    plan_path='metadata/baliphy_sampler_resource_observation_plan_20261003.json'
    output='results/ancestral/full-baliphy-sampler-resource-observation-20261003-v1'
    plan=dict(sampler_plan=str(native_path),sampler_launch=original_path,output=output,pins=pins,
        resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=1,memory_gib=2,swap_gib=0,
            reader_cpus=2,reader_memory_gib=4,poll_seconds=1,blas_threads=1,minimum_free_disk_gib=64,
            output_allowance_gib=16,runtime_depends_on_original_controller=True,runtime_uncalibrated=True,finish_eta=None,
            planning_observation_days=[.1,4],native_jobs_started=0,gpu=False,new_cost_usd=0,
            available_memory_gib=psutil.virtual_memory().available/2**30,
            available_disk_gib=shutil.disk_usage('.').free/2**30,
            caveat='Continuous cadence includes processing/scheduling delays; actual gaps retained. Streaming reader avoids holding the event history in memory. Allowance is a planning estimate, not a hard output quota or calibrated run ETA.'),
        completion='metadata/baliphy_sampler_resource_observation_completed_20261003.json',
        launch_inventory='metadata/baliphy_sampler_resource_observation_launches_20261003.json',scope=scope)
    create(plan_path,plan)
    producer=launch('baliphy-sampler-resource-observation',
        [sys.executable,'scripts/run_baliphy_sampler_resource_observation.py','--plan',plan_path],[],plan_path,cpus=1,memory=2)
    reader=launch('baliphy-sampler-resource-observation-readback',
        [sys.executable,'scripts/run_baliphy_sampler_resource_observation.py','--plan',plan_path,'--reader'],[producer],plan_path,cpus=2,memory=4)
    completion_path='metadata/baliphy_sampler_resource_observation_completion_plan_20261003.json'
    spec=dict(source_plan=plan_path,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,completed_status=COMPLETED_STATUS,
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={p:sha(p) for p in [plan_path,producer,reader,'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=scope)
    create(completion_path,spec)
    closer=launch('baliphy-sampler-resource-observation-closure',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',completion_path],
        [producer,reader],completion_path,cpus=2,memory=4)
    create(plan['launch_inventory'],dict(status='launched_full_read_only_original_sampler_resource_observer',
        source_plan=plan_path,source_plan_sha256=sha(plan_path),launches=[producer,reader,closer],
        original_sampler_launch=original_path,new_native_jobs=0,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
