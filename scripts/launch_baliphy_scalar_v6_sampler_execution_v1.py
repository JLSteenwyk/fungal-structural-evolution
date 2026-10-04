#!/usr/bin/env python3
"""Queue complete V6 computational sampling and concurrent read-only telemetry."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_scalar_v6_sampler_v1 import PRODUCER_STATUS, READER_STATUS, COMPLETED_STATUS, SUMMARY_FIELDS
from run_baliphy_scalar_v6_resource_observer_v1 import qualified_jobs, SUMMARY_FIELDS as OBSERVER_FIELDS
import run_baliphy_scalar_v6_resource_observer_v1 as observer


def create(path, value):
    with Path(path).open('x') as f: json.dump(value, f, indent=2); f.write('\n')


def launch(label, command, dependencies, source, cpus=2, memory=32):
    prefix = 'metadata/baliphy_scalar_v6_execution_' + label + '_20261004_v1'
    unit = 'fungal-scalar-v6-execution-' + label + '-20261004-v1.service'
    wait_path = prefix + '_wait_plan.json'; launch_path = prefix + '_launch.json'
    assert not Path(launch_path).exists()
    bindings = {path: sha(path) for path in [source, *dependencies, str(Path(__file__)),
        'scripts/run_after_verified_dependencies_v2.py']}
    create(wait_path, dict(dependencies=dependencies, command=command, pins=bindings))
    cmd = ['systemd-run', '--user', '--collect', '--unit=' + unit,
        '--working-directory=' + str(Path.cwd()), '-p', 'CPUQuota='+str(cpus*100)+'%',
        '-p', 'MemoryMax='+str(memory)+'G', '-p', 'MemorySwapMax=0', '-p', 'TasksMax=512',
        '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1',
        '--setenv=OMP_NUM_THREADS=1', '--setenv=MKL_NUM_THREADS=1',
        sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    subprocess.run(cmd, check=True)
    pid = 0
    for _ in range(20):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', unit,
            '-p', 'MainPID', '--value'], text=True))
        if pid: break
        time.sleep(.1)
    assert pid > 0
    proc = psutil.Process(pid)
    expected = [sys.executable, 'scripts/run_after_verified_dependencies_v2.py', '--plan', wait_path]
    assert proc.cmdline() == expected
    cg = subprocess.check_output(['systemctl', '--user', 'show', unit,
        '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {k: (Path('/sys/fs/cgroup') / cg.lstrip('/') / k).read_text().strip()
        for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': str(cpus*100000)+' 100000', 'memory.max': str(memory * 2**30), 'memory.swap.max': '0'}
    create(launch_path, dict(unit=unit, pid=pid, created=proc.create_time(), cmdline=expected,
        plan=wait_path, plan_sha256=sha(wait_path), checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits, systemd_launch_command=cmd))
    print(json.dumps(dict(launch=launch_path, pid=pid, unit=unit)), flush=True)
    return launch_path



def validate(plan_path):
    plan=json.loads(Path(plan_path).read_text());verify(plan['pins'])
    assert plan['scalar_schema']=='native-cjson-explicit-special-values-v6'
    jobs=qualified_jobs(plan);assert len(jobs)==1620
    for field, expected in [('sampler_software_transport','verified_original_scalar_v6_sampler_controller_software_wait_zero'),
                            ('observer_software_transport','verified_original_scalar_v6_resource_observer_software_wait_zero')]:
        proof=json.loads(Path(plan[field]).read_text())
        assert proof['status']==expected and proof['actual_tool_terminal_exit_code']==0
        assert proof['entire_terminal_payload_matched'];verify(proof['source_hashes'])
    r=plan['resources']
    assert (r['workers'],r['cpus'],r['memory_gib'],r['reservation_capacity_gib'],r['swap_gib'],r['blas_threads'],r['iterations'])==(16,16,200,192,0,1,20)
    assert r['posterior_qualified'] is False
    assert len(plan['dependencies'])==5 and plan['dependencies'][0]==plan['startup_closure_launch']
    previous=json.loads(Path(plan['previous_joint_sampler_plan']).read_text())
    assert plan['dependencies'][1:]==previous['dependencies']
    assert psutil.virtual_memory().available >= r['memory_gib']*2**30
    assert psutil.disk_usage('.').free >= r['minimum_free_disk_gib']*2**30
    assert not Path(plan['output']).exists() and not Path(plan['observer_output']).exists()
    for path in plan['dependencies']:
        record=json.loads(Path(path).read_text());record['launch']=path
        assert sha(record['plan'])==record['plan_sha256']
        if fingerprint(record) is None:journal_terminal(record)
    return plan


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--validate-only',action='store_true')
    args=parser.parse_args();plan=validate(args.plan)
    if args.validate_only:
        print('Verified full V6 software/jobs/resources and exact prerequisite identities; no launch.')
        return
    producer=launch('sampler-producer',[sys.executable,'scripts/run_baliphy_scalar_v6_sampler_v1.py','--plan',str(args.plan)],
        plan['dependencies'],str(args.plan),cpus=16,memory=200)
    reader=launch('sampler-reader',[sys.executable,'scripts/run_baliphy_scalar_v6_sampler_v1.py','--plan',str(args.plan),'--reader'],
        [producer],str(args.plan))
    root=Path(plan['output'])
    completion=dict(source_plan=str(args.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,completed_status=COMPLETED_STATUS,
        summary_fields=SUMMARY_FIELDS+['reservation_audit'],launches=[producer,reader],
        pins={p:sha(p) for p in [str(args.plan),producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=plan['scope'])
    create(plan['completion_plan'],completion)
    closer=launch('sampler-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['completion_plan']],
        [producer,reader],plan['completion_plan'])
    observer_plan=dict(sampler_plan=str(args.plan),sampler_launch=producer,output=plan['observer_output'],
        resources=dict(cpus=1,memory_gib=2,reader_cpus=2,reader_memory_gib=4,swap_gib=0,poll_seconds=1,minimum_free_disk_gib=64,
            output_allowance_gib=16,finish_eta=None,gpu=False,new_cost_usd=0),
        completion=plan['observer_completion'],pins={p:sha(p) for p in [str(args.plan),producer,
            'scripts/run_baliphy_scalar_v6_resource_observer_v1.py','scripts/baliphy_sampler_resource_observation.py',
            plan['observer_software_transport'],str(Path(__file__))]},
        scope='Full1620V6role read-only exact original sampler/native telemetry. Missing/fast/transitioning/failure attempts explicit; no exact per-native peak guarantee. Source/artifact/independent serialized reader/two original observer journals close after sampler termination. No retries, native controls, posterior or biological acceptance.')
    create(plan['observer_plan'],observer_plan)
    # Observe the producer while it is still waiting on prerequisite closure;
    # dependency is startup only, never producer completion. Fast misses remain explicit.
    telemetry=launch('observer-producer',[sys.executable,'scripts/run_baliphy_scalar_v6_resource_observer_v1.py','--plan',plan['observer_plan']],
        [plan['startup_closure_launch']],plan['observer_plan'],cpus=1,memory=2)
    telemetry_reader=launch('observer-reader',[sys.executable,'scripts/run_baliphy_scalar_v6_resource_observer_v1.py','--plan',plan['observer_plan'],'--reader'],
        [telemetry],plan['observer_plan'],cpus=2,memory=4)
    obroot=Path(plan['observer_output'])
    obclosure=dict(source_plan=plan['observer_plan'],producer_receipt=str(obroot/'receipt.json'),independent_readback=str(obroot/'readback.json'),
        producer_status=observer.PRODUCER_STATUS,reader_status=observer.READER_STATUS,completed_status=observer.COMPLETED_STATUS,
        summary_fields=OBSERVER_FIELDS,launches=[telemetry,telemetry_reader],
        pins={p:sha(p) for p in [plan['observer_plan'],telemetry,telemetry_reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['observer_completion'],scope=observer_plan['scope'])
    create(plan['observer_completion_plan'],obclosure)
    telemetry_closer=launch('observer-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['observer_completion_plan']],
        [telemetry,telemetry_reader],plan['observer_completion_plan'],cpus=2,memory=4)
    create(plan['launch_inventory'],dict(status='queued_full_scalar_v6_short_sampler_and_read_only_observer',
        source_plan=str(args.plan),source_plan_sha256=sha(args.plan),sampler_launches=[producer,reader,closer],
        observer_launches=[telemetry,telemetry_reader,telemetry_closer],original_prerequisite_closures=plan['dependencies'],
        expected_native_roles=1620,iterations_per_role=20,posterior_qualified=False,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
