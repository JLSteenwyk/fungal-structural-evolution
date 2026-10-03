#!/usr/bin/env python3
"""Exercise full telemetry serialization with explicitly artificial observations."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

from ancestral_chain_attempt import sha, write_json
from baliphy_sampler_resource_observation import describe, probe_native, SUMMARY_FIELDS
from record_project_runtime_checkpoint_v4 import fingerprint
from record_baliphy_sampler_resource_observation_checkpoint import last_event
import resource
import time
import run_baliphy_joint_sampler_resource_observation as workflow
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError):return
    raise AssertionError('Altered telemetry export or completed restart accepted')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--native-validation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();native_gate=json.loads(args.native_validation.read_text());verify(dict(pins=native_gate['source_hashes']))
    started=time.monotonic()
    root=args.output.resolve();root.mkdir(exist_ok=False)
    original_path=Path('metadata/baliphy_joint_sampler_qualification_v2_plan_20261003.json')
    original=json.loads(original_path.read_text());verify(original)
    jobs=json.loads(Path(original['jobs']).read_text());assert len(jobs)==1620
    for job in jobs:
        assert job['config']['command'][job['config']['command'].index('--seed')+1]==str(job['chain']['seed'])
        assert job['memory_reservation_bytes'] in [12*2**30,48*2**30]
        for path,h in job['config']['pins'].items():assert sha(path)==h,path
    native_root=root/'artificial-native';native_root.mkdir();native_plan=root/'artificial-native-plan.json'
    write_json(native_plan,dict(output=str(native_root),jobs=original['jobs'],pins={}))
    launch_path=root/'artificial-launch.json';record=dict(pid=321,created=100.,cmdline=['artificial-controller'],
        actual_cgroup_limits={'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'})
    write_json(launch_path,record);terminal=dict(launch=str(launch_path),pid=321,created=100.,
        cmdline=record['cmdline'],status='verified_original_terminal_success_with_bound_completed_artifacts',
        scope='Software fixture only: original invocation and terminal journal are mocked')
    plan_path=root/'plan.json';plan=dict(output=str(root/'stage'),sampler_plan=str(native_plan),
        sampler_launch=str(launch_path),resources=dict(poll_seconds=1,minimum_free_disk_gib=1),pins={},
        scope='Artificial full telemetry serialization fixture; native execution and journal proof mocked')
    write_json(plan_path,plan)
    def fake_capture(plan,jobs,stage):
        descriptions=[];samples=[]
        for i,j in enumerate(jobs):
            cid=j['chain']['chain_id'];folder=native_root/'attempts'/cid/'attempt-0001';folder.mkdir(parents=True)
            write_json(folder.parent/'configuration.json',j['config']);write_json(folder/'command.json',j['config']['command'])
            write_json(folder/'process.json',dict(pid=1000+i,pgid=1000+i,created=100.,command=j['config']['command']))
            descriptor=describe(j,folder/'process.json');descriptions.append(descriptor)
            digest=hashlib.sha256(json.dumps(j['config'],sort_keys=True,separators=(',', ':'),allow_nan=False).encode()).hexdigest()
            write_json(folder/'receipt.json',dict(configuration_sha256=digest,exit_code=1 if i in [0,4] else 0,
                status='failed' if i in [0,4] else 'exited_zero_pending_scientific_validation',
                artifacts={name:sha(folder/name) for name in ['command.json','process.json']}))
            if i<1618:
                samples.append(dict(chain_id=cid,pid=1000+i,created=100.,identity_sha256=descriptor['identity_sha256'],
                    checked_utc='2026-10-03T00:00:00+00:00',status='verified_live_native_resource_observation',
                    cgroup='/artificial-unit',actual_limits={'Max address space':j['memory_reservation_bytes'],
                        'Max cpu time':j['config']['timeout_seconds'],'Max file size':2*2**30},
                    reported_memory_bytes={'VmPeak':2000,'VmSize':1000,'VmHWM':800,'VmRSS':600,'VmSwap':0},cpu_seconds=3.))
        snapshot=dict(sequence=0,sampler_controller_pid=321,sampler_controller_created=100.,
            elapsed_since_previous_snapshot_seconds=1.,new_attempts=descriptions,native_observations=samples,
            cgroup=dict(path='/artificial-unit',limits=record['actual_cgroup_limits'],
                memory_bytes={'memory.current':10000,'memory.peak':20000}))
        (stage/'observations.jsonl').write_text(json.dumps(snapshot)+'\n')
        write_json(stage/'attempts_at_termination.json',descriptions);write_json(stage/'sampler_terminal.json',terminal)
        write_json(native_root/'receipt.json',dict(full_chains=1620,plan_sha256=sha(native_plan),scientific_eligibility=False))
    negatives=[]
    with patch.object(workflow,'runtime_caps',return_value={'artificial_fixture':True}),patch.object(
        workflow,'capture',side_effect=fake_capture),patch.object(workflow,'fingerprint',return_value=None),patch.object(
        workflow,'journal_terminal',return_value=terminal):
        producer=workflow.run(plan_path)
        assert producer['roles_with_native_attempt']==1620 and producer['roles_with_live_observation']==1618
        assert producer['roles_without_live_observation']==2 and producer['native_outcome_status_counts']['failed']==2
        rejected(lambda:workflow.run(plan_path));negatives.append('completed_producer_restart')
        exported=Path(plan['output'])/'roles.json';original_bytes=exported.read_bytes();rows=json.loads(original_bytes)
        for name in ['missing_role','duplicated_role','changed_alias','changed_observation_count','invented_missing_peak','scientific_acceptance']:
            bad=copy.deepcopy(rows)
            if name=='missing_role':bad.pop()
            elif name=='duplicated_role':bad[1]=copy.deepcopy(bad[0])
            elif name=='changed_alias':bad[1]['original_configuration_ids']=['invented']
            elif name=='changed_observation_count':bad[1]['live_observations']+=1
            elif name=='invented_missing_peak':bad[-1]['maximum_reported_memory_bytes']['VmHWM']=99999999
            else:bad[1]['scientific_eligibility']=True
            write_json(exported,bad);rejected(lambda:workflow.run(plan_path,reader=True));negatives.append(name)
            exported.write_bytes(original_bytes)
        reader=workflow.run(plan_path,reader=True);assert all(producer[k]==reader[k] for k in SUMMARY_FIELDS)
        rejected(lambda:workflow.run(plan_path,reader=True));negatives.append('completed_reader_restart')
    # Reuse the qualified shared procfs decoder, then passively inspect current
    # historical native workers. These are not new joint-stage measurements.
    historical_path=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    historical=json.loads(historical_path.read_text());verify(historical)
    historical_jobs={job['chain']['chain_id']:job for job in json.loads(Path(historical['jobs']).read_text())}
    observed=[];descriptors=[]
    historical_journal=Path('results/ancestral/full-baliphy-sampler-resource-observation-20261003-v1/observations.jsonl')
    for _ in range(20):
        latest=last_event(historical_journal)
        for entry in latest['native_observations']:
            if entry['status']!='verified_live_native_resource_observation':continue
            job=historical_jobs[entry['chain_id']]
            identity=Path(historical['output'])/'attempts'/entry['chain_id']/'attempt-0001'/'process.json'
            descriptor=describe(job,identity)
            sample=probe_native(job,descriptor,latest['cgroup']['path'])
            descriptors.append(descriptor);observed.append(sample)
        if any(sample['status']=='verified_live_native_resource_observation' for sample in observed):break
        time.sleep(.1)
    live=[sample for sample in observed if sample['status']=='verified_live_native_resource_observation']
    assert live,'Original native workers must be revalidated, not assumed live'
    live_path=root/'historical_live_native_rechecks.json';write_json(live_path,dict(descriptors=descriptors,observations=observed))
    joint_launch=Path('metadata/baliphy_joint_sampler_qualification_v2_launch_20261003.json')
    joint=json.loads(joint_launch.read_text());assert fingerprint(joint) is not None
    assert joint['actual_cgroup_limits']=={'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'}
    usage=resource.getrusage(resource.RUSAGE_SELF)
    result=dict(native_gate);result.update(status='passed_full_joint_sampler_resource_observation_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        historical_qualified_native_probe_readings=native_gate['actual_live_native_observations'],
        actual_live_native_observations=len(live),actual_live_observation_source='Original historical sampler only; joint sampler remains queued',
        original_joint_controller_pid=joint['pid'],full_joint_configurations_checked=1620,new_native_inference_runs=0,
        elapsed_seconds=time.monotonic()-started,self_cpu_seconds=usage.ru_utime+usage.ru_stime,
        self_peak_rss_bytes=usage.ru_maxrss*1024,
        scope='All1620jointconfigs, shared immutable procfs decoder requalified by prior native fixture gate, current original historical native workers passively probed; full artificial joint-stage serialization/missing/failure accounting with mocked capture and original journals. New joint controller verified live/queued. Not actual joint native observations, production closure, precise final native peaks, long-chain resources or posterior acceptance.',

        native_validation=str(args.native_validation),native_validation_sha256=sha(args.native_validation),
        full_producer_and_reader_serialization_checked=True,artificial_failed_roles_retained=2,
        artificial_missing_live_observations_retained=2,serialization_rejection_cases=negatives,
        serialization_scope='Full1620real configurations with artificial native attempts/receipts/procfs samples and mocked original-controller journals. Real full producer/reader/replay/hash/export code exercised. This is software accounting proof, not actual production measurements or runtime closure.')
    result['source_hashes']=dict(native_gate['source_hashes'])
    for p in [Path(__file__),args.native_validation,original_path,Path(original['jobs']),plan_path,native_plan,launch_path,
              historical_path,Path(historical['jobs']),live_path,joint_launch,Path('scripts/run_baliphy_joint_sampler_resource_observation.py')]:result['source_hashes'][str(p)]=sha(p)
    for descriptor in descriptors:
        for field in ['identity','configuration','command']:
            path=descriptor[field+'_path'];result['source_hashes'][path]=descriptor[field+'_sha256']
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],full_roles_checked=1620,
        full_producer_and_reader_serialization_checked=True,rejection_cases=negatives)),flush=True)


if __name__=='__main__':main()
