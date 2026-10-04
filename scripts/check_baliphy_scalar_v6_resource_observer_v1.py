#!/usr/bin/env python3
"""Qualify V6 full-grid resource observation with explicit artificial telemetry."""
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
import run_baliphy_scalar_v6_resource_observer_v1 as workflow
from baliphy_scalar_json_logger_v6c import SCHEMA
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind
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
    original_path=Path('metadata/baliphy_scalar_v6_sampler_controller_software_validation_20261004_v1.json')
    original=json.loads(original_path.read_text());verify(dict(pins=original['source_hashes']))
    job_path=Path('data/software_audits/baliphy-scalar-json-v6-sampler-controller-20261004-v1/actual_full_grid_jobs.json').resolve()
    jobs=workflow.qualified_jobs(dict(jobs=str(job_path),scalar_schema=SCHEMA));assert len(jobs)==1620
    for job in jobs:
        assert job['config']['command'][job['config']['command'].index('--seed')+1]==str(job['chain']['seed'])
        assert job['memory_reservation_bytes'] in [12*2**30,48*2**30]
        for path,h in job['config']['pins'].items():assert sha(path)==h,path
    native_root=root/'artificial-native';native_root.mkdir();native_plan=root/'artificial-native-plan.json'
    write_json(native_plan,dict(output=str(native_root),jobs=str(job_path),scalar_schema=SCHEMA,pins={}))
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
        write_json(native_root/'receipt.json',dict(status='complete_full_scalar_v6_short_sampler_dispositions_pending_readback_v1',full_chains=1620,plan_sha256=sha(native_plan),scientific_eligibility=False))
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
    # Reuse the exact hash-bound native procfs qualification: the original
    # native probe is retained evidence, not a claim of current native work.
    assert native_gate['actual_live_native_observations'] > 0
    assert sha('scripts/baliphy_sampler_resource_observation.py') == native_gate['source_hashes']['scripts/baliphy_sampler_resource_observation.py']
    assert native_gate['full_missing_observation_accounting_checked']
    assert native_gate['empty_or_transitioning_argv_not_accepted_as_native']
    assert native_gate['stable_unrelated_command_rejected'] and native_gate['wrapper_not_mistaken_for_native']
    descriptors = []; observed = []; live = []
    # Admission rejects foreign schema and any changed role, including aliases
    # and prior metadata, before telemetry output is created. Private copies only.
    design_negatives=[]
    for case in ['foreign_schema','changed_seed','changed_prior','changed_alias','missing_role','duplicate_role']:
        altered=copy.deepcopy(jobs); native=dict(jobs=str(job_path),scalar_schema=SCHEMA)
        if case=='foreign_schema':native['scalar_schema']='old-unqualified-schema'
        else:
            if case=='changed_seed':altered[0]['chain']['seed']+=1
            elif case=='changed_prior':altered[0]['chain']['prior_label']='foreign'
            elif case=='changed_alias':altered[0]['chain']['original_configuration_ids']=['foreign']
            elif case=='missing_role':altered.pop()
            else:altered[1]=copy.deepcopy(altered[0])
            altered_path=root/(case+'-jobs.json');write_json(altered_path,altered);native['jobs']=str(altered_path)
        rejected(lambda:workflow.qualified_jobs(native));design_negatives.append(case)
    usage=resource.getrusage(resource.RUSAGE_SELF)
    result=dict(native_gate);result.update(status='passed_full_scalar_v6_sampler_resource_observer_software_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        historical_qualified_native_probe_readings=native_gate['actual_live_native_observations'],
        actual_live_native_observations=len(live),actual_live_observation_source='No current native readings claimed; exact preserved shared decoder native qualification reused',
        current_read_only_controller_verified=False,future_native_observer_launched=False,full_v6_configurations_checked=1620,design_rejection_cases=design_negatives,new_native_inference_runs=0,
        elapsed_seconds=time.monotonic()-started,self_cpu_seconds=usage.ru_utime+usage.ru_stime,
        self_peak_rss_bytes=usage.ru_maxrss*1024,
        scope='All1620jointconfigs, shared immutable procfs decoder requalified by prior native fixture gate, preserved original native procfs fixture qualification hash-bound and reused; zero current native readings claimed; full artificial joint-stage serialization/missing/failure accounting with mocked capture and original journals. No future native observer/controller claimed live; exact future V6 job bytes bound. Not actual joint native observations, production closure, precise final native peaks, long-chain resources or posterior acceptance.',

        native_validation=str(args.native_validation),native_validation_sha256=sha(args.native_validation),
        full_producer_and_reader_serialization_checked=True,artificial_failed_roles_retained=2,
        artificial_missing_live_observations_retained=2,serialization_rejection_cases=negatives,
        serialization_scope='Full1620real configurations with artificial native attempts/receipts/procfs samples and mocked original-controller journals. Real full producer/reader/replay/hash/export code exercised. This is software accounting proof, not actual production measurements or runtime closure.')
    result['source_hashes']=dict(native_gate['source_hashes'])
    for p in [Path(__file__),args.native_validation,original_path,job_path,plan_path,native_plan,launch_path,
              Path('scripts/run_baliphy_scalar_v6_resource_observer_v1.py')]:result['source_hashes'][str(p)]=sha(p)
    for descriptor in descriptors:
        for field in ['identity','configuration','command']:
            path=descriptor[field+'_path'];result['source_hashes'][path]=descriptor[field+'_sha256']
    for n,h in original['source_hashes'].items():bind(result['source_hashes'],n,h)
    modules=project_sources(result['source_hashes'],[Path(__file__),Path('scripts/run_baliphy_scalar_v6_resource_observer_v1.py')])
    result['transitive_project_source_modules']=len(modules)
    for path in root.rglob('*'):
        if path.is_file():bind(result['source_hashes'],path)
    verify(dict(pins=result['source_hashes']))
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],full_roles_checked=1620,
        full_producer_and_reader_serialization_checked=True,rejection_cases=negatives)),flush=True)


if __name__=='__main__':main()
