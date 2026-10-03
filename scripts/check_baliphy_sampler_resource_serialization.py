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
from baliphy_sampler_resource_observation import describe, SUMMARY_FIELDS
import run_baliphy_sampler_resource_observation as workflow
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError):return
    raise AssertionError('Altered telemetry export or completed restart accepted')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--native-validation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();native_gate=json.loads(args.native_validation.read_text());verify(dict(pins=native_gate['source_hashes']))
    root=args.output.resolve();root.mkdir(exist_ok=False)
    original_path=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    original=json.loads(original_path.read_text());jobs=json.loads(Path(original['jobs']).read_text());assert len(jobs)==1620
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
    result=dict(native_gate);result.update(checked_utc=datetime.now(timezone.utc).isoformat(),
        native_validation=str(args.native_validation),native_validation_sha256=sha(args.native_validation),
        full_producer_and_reader_serialization_checked=True,artificial_failed_roles_retained=2,
        artificial_missing_live_observations_retained=2,serialization_rejection_cases=negatives,
        serialization_scope='Full1620real configurations with artificial native attempts/receipts/procfs samples and mocked original-controller journals. Real full producer/reader/replay/hash/export code exercised. This is software accounting proof, not actual production measurements or runtime closure.')
    result['source_hashes']=dict(native_gate['source_hashes'])
    for p in [Path(__file__),args.native_validation,original_path,plan_path,native_plan,launch_path]:result['source_hashes'][str(p)]=sha(p)
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],full_roles_checked=1620,
        full_producer_and_reader_serialization_checked=True,rejection_cases=negatives)),flush=True)


if __name__=='__main__':main()
