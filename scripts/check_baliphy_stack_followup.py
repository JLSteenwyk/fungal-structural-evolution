#!/usr/bin/env python3
"""Exercise full-role accounting with explicit native and source-gate mocks."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import time
from unittest.mock import patch

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_sampler_qualification_v3 import SUCCESS
from baliphy_sampler_resource_observation import describe
from baliphy_stack_followup import build_jobs, forbidden_seeds, SUMMARY_FIELDS
from reference_measurement_union_sources import verify
import run_baliphy_stack_followup as workflow


def rejected(action):
    try:action()
    except (AssertionError,ValueError,KeyError,RuntimeError,FileNotFoundError):return
    raise AssertionError('Altered follow-up accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    root=a.output.resolve();root.mkdir(exist_ok=False);assert not a.receipt.exists();started=time.monotonic()
    sp=Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    source=json.loads(sp.read_text());originals=json.loads(Path(source['jobs']).read_text())
    rows_path=Path(source['output'])/'dispositions.json';rows=json.loads(rows_path.read_text())
    hp=Path(source['historical_sampler_plan']);history=json.loads(hp.read_text());historical=json.loads(Path(history['jobs']).read_text())
    dp=Path('metadata/baliphy_stack_followup_diagnosis_20261003_v1.json');diagnosis=json.loads(dp.read_text())
    verify(diagnosis['source_hashes']);jobs=build_jobs(originals,rows,historical,diagnosis)
    assert len(jobs)==24 and {j['chain']['seed'] for j in jobs}.isdisjoint(forbidden_seeds(originals,historical))
    design_negatives=[]
    for name in ['missing_original_role','duplicate_original_role','changed_original_alias','changed_original_prior',
                 'unexpected_failure_family','different_signal','allocation_failure','nonempty_failed_frame','missing_diagnosed_input']:
        oo=copy.deepcopy(originals);rr=copy.deepcopy(rows);dd=copy.deepcopy(diagnosis)
        failure=next(r for r in rr if r['exit_code']!=0)
        if name=='missing_original_role':oo.pop()
        elif name=='duplicate_original_role':oo[1]=copy.deepcopy(oo[0])
        elif name=='changed_original_alias':next(j for j in oo if j['chain']['chain_id']==failure['chain_id'])['chain']['original_configuration_ids']=[]
        elif name=='changed_original_prior':oo[0]['chain']['prior_label']='invented'
        elif name=='unexpected_failure_family':next(j for j in oo if j['chain']['chain_id']==failure['chain_id'])['chain']['family']='other'
        elif name=='different_signal':failure['exit_code']=-9
        elif name=='allocation_failure':failure['bad_alloc']=True
        elif name=='nonempty_failed_frame':failure['joint_frames']=[{'invented':True}]
        else:dd['failed_input_group_counts'].pop(next(iter(dd['failed_input_group_counts'])))
        rejected(lambda:build_jobs(oo,rr,historical,dd));design_negatives.append(name)
    negatives=[]
    rejected(lambda:workflow.stack_limit('Max stack size 8388608 8388608 bytes\n'));negatives.append('old_stack_limit')
    rejected(lambda:workflow.stack_limit('Max stack size 67108864 8388608 bytes\n'));negatives.append('different_hard_stack_limit')
    assert workflow.stack_limit('Max stack size 67108864 67108864 bytes\n')==64*2**20
    jp=root/'jobs.json';write_json(jp,jobs);pp=root/'plan.json';stage=root/'pipeline'
    plan=dict(source_plan=str(sp),source_resource_plan='software-mocked',diagnosis=str(dp),jobs=str(jp),
        mapping=source['mapping'],output=str(stage),pins={str(jp):sha(jp)},
        resources=dict(cpus=4,memory_gib=200,reservation_capacity_gib=192,workers=4,minimum_free_disk_gib=1,
                       iterations=20,native_stack_bytes=64*2**20),scope='Explicit synthetic native attempts and mocked source closure/cgroup caps; no scientific admission.')
    write_json(pp,plan);fake_rows={};first=next(r for r in rows if r['status']==SUCCESS)
    artificial_failures={jobs[i]['chain']['chain_id'] for i in [0,4,8]}
    def fake_execute(job,root,digest,mapping,budget,minimum_disk):
        c=job['chain'];cid=c['chain_id']
        with budget.reserve(cid,job['memory_reservation_bytes']):
            folder=root/'attempts'/cid/'attempt-0001';folder.mkdir(parents=True)
            write_json(folder.parent/'configuration.json',job['config']);write_json(folder/'command.json',job['config']['command'])
            write_json(folder/'process.json',dict(pid=10000+jobs.index(job),pgid=10000+jobs.index(job),created=100.,command=job['config']['command']))
            encoded=json.dumps(job['config'],sort_keys=True,separators=(',',':'),allow_nan=False)
            failed=cid in artificial_failures
            rp=folder/'receipt.json';write_json(rp,dict(configuration_sha256=hashlib.sha256(encoded.encode()).hexdigest(),
                exit_code=1 if failed else 0,status='failed' if failed else 'exited_zero_pending_scientific_validation',
                artifacts={p.name:sha(p) for p in [folder/'process.json',folder/'command.json']}))
            row=copy.deepcopy(first)
            row.update(chain_id=cid,effective_input_group=c['effective_input_group'],model_input_identity=c['effective_input_group']+'-'+c['prior_label'],
                prior_label=c['prior_label'],chain_role=c['chain'],family=c['family'],original_configuration_ids=c['original_configuration_ids'],
                seed=c['seed'],source_seed=job['source_seed'],native_receipt=str(rp),native_receipt_sha256=sha(rp),plan_sha256=digest,
                exit_code=1 if failed else 0,native_status='failed' if failed else 'exited_zero_pending_scientific_validation')
            if failed:
                row.update(status='unsuccessful_sampler_qualification_attempt_retained',saved_alignments=0,candidate_frames=0,
                           joint_frames=[],ancestral_categories_available=False,same_record_sequence_category_correspondence=False)
            else:
                export=root/'frames'/cid;export.mkdir(parents=True)
                for frame in row['joint_frames']:
                    destination=export/Path(frame['projection_array']).name
                    shutil.copyfile(frame['projection_array'],destination)
                    frame['projection_array']=str(destination);frame['projection_array_sha256']=sha(destination)
            fake_rows[cid]=row;write_json(root/'chains'/(cid+'.json'),row)
            return row
    def fake_capture(root,jobs,controller,stop):
        stop.wait()
        samples=[]
        for j in jobs:
            cid=j['chain']['chain_id'];d=describe(j,root/'attempts'/cid/'attempt-0001/process.json')
            samples.append(dict(chain_id=cid,pid=d['pid'],created=d['created'],identity_sha256=d['identity_sha256'],
                status='verified_live_native_resource_observation',command=j['config']['command'][j['config']['command'].index('--')+1:],
                checked_utc='2026-10-03T00:00:00+00:00',native_stack_bytes=64*2**20,cgroup='/synthetic',
                actual_limits={'Max address space':48*2**30,'Max cpu time':j['config']['timeout_seconds'],'Max file size':2*2**30},
                reported_memory_bytes={'VmPeak':2000,'VmSize':1000,'VmHWM':800,'VmRSS':600,'VmSwap':0},cpu_seconds=1.))
        event=dict(sequence=0,controller_pid=controller['pid'],controller_created=controller['created'],
            elapsed_since_previous_snapshot_seconds=1.,cgroup=dict(path='/synthetic',limits=controller['actual_cgroup_limits'],
            memory_bytes={'memory.current':100000,'memory.peak':200000,'memory.swap.current':0},
            memory_events={'low':0,'high':0,'max':0,'oom':0,'oom_kill':0,'oom_group_kill':0}),native_observations=samples)
        (root/'observations.jsonl').write_text(json.dumps(event)+'\n')
    def fake_inspect(job,rp,digest,mapping,export,allow_export_creation=False):
        assert allow_export_creation is False
        row=fake_rows[job['chain']['chain_id']]
        for frame in row['joint_frames']:assert sha(frame['projection_array'])==frame['projection_array_sha256']
        return copy.deepcopy(row)
    def fake_caps(plan,reader=False):
        return {'cpu.max':'200000 100000' if reader else '400000 100000',
                'memory.max':str((32 if reader else 200)*2**30),'memory.swap.max':'0'}
    with patch.object(workflow,'prerequisites',return_value=(originals,rows,jobs,{})),patch.object(
        workflow,'runtime_caps',side_effect=fake_caps),patch.object(workflow,'execute_job',side_effect=fake_execute),patch.object(
        workflow,'capture',side_effect=fake_capture),patch.object(workflow,'inspect',side_effect=fake_inspect),patch.dict(
        'os.environ',{'INVOCATION_ID':'explicit-synthetic-software-invocation'}):
        producer=workflow.run(pp)
        assert producer['followup_roles']==24 and producer['original_failed_roles']==24
        assert producer['selected_checked_sampler_attempts']==1617 and producer['selected_unsuccessful_sampler_attempts']==3
        assert producer['reservation_audit']['ledger_events']==48
        assert producer['resource_audit']['roles_with_live_stack_limit_observation']==24
        rejected(lambda:workflow.run(pp));negatives.append('completed_producer_restart')
        mutable=stage/'full_role_accounting.json';saved=mutable.read_bytes();altered=json.loads(saved);altered.pop()
        mutable.write_text(json.dumps(altered));rejected(lambda:workflow.run(pp,True));negatives.append('missing_original_accounting_role');mutable.write_bytes(saved)
        array=next(iter((stage/'frames').glob('*/*')));saved_array=array.read_bytes();array.unlink()
        rejected(lambda:workflow.run(pp,True));assert not array.exists();negatives.append('missing_export_not_recreated');array.write_bytes(saved_array)
        obs=stage/'observations.jsonl';saved_obs=obs.read_bytes()
        for name in ['wrong_live_stack','changed_live_command','cgroup_oom_event','changed_observation_identity']:
            altered=json.loads(saved_obs)
            if name=='wrong_live_stack':altered['native_observations'][0]['native_stack_bytes']=8*2**20
            elif name=='changed_live_command':altered['native_observations'][0]['command']=['invented']
            elif name=='cgroup_oom_event':altered['cgroup']['memory_events']['oom']=1
            else:altered['native_observations'][0]['created']=101.
            obs.write_text(json.dumps(altered)+'\n')
            rejected(lambda:workflow.resource_replay(stage,jobs,fake_caps(plan)));negatives.append(name)
        obs.write_bytes(saved_obs)
        reader=workflow.run(pp,True);assert all(reader[k]==producer[k] for k in SUMMARY_FIELDS)
        rejected(lambda:workflow.run(pp,True));negatives.append('completed_reader_restart')
    bindings={str(p):sha(p) for p in [Path(__file__),sp,Path(source['jobs']),rows_path,hp,Path(history['jobs']),dp,
        Path('scripts/baliphy_stack_followup.py'),Path('scripts/run_baliphy_stack_followup.py'),pp,jp]}
    for p in stage.rglob('*'):
        if p.is_file():bindings[str(p)]=sha(p)
    verify(bindings)
    result=dict(status='passed_full_role_stack_followup_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_original_roles=1620,full_original_quartets=405,full_original_inputs=135,full_original_aliases=324,
        followup_roles=24,followup_priors=3,followup_chain_roles=4,diagnosed_inputs=2,
        fresh_disjoint_seeds=24,design_cases_rejected=design_negatives,serialization_and_resource_cases_rejected=negatives,
        artificial_followup_failures_retained=3,artificial_selected_successful_roles=1617,
        full_producer_reader_serialization_passed=True,native_execution_mocked=True,source_closure_mocked=True,
        cgroup_and_procfs_observations_mocked=True,scientific_eligibility=False,posterior_qualified=False,
        elapsed_seconds=time.monotonic()-started,source_hashes=bindings,
        scope='Complete1620real-role metadata and24generated follow-up jobs through full serialization with explicit synthetic native/telemetry/source-gate fixtures. Copied original arrays are software fixtures only. All historical failures retained; malformed accounting/exports/limits/identity/OOM and restart rejected. This proves software contracts, not24native stack-corrected sampling successes or actual new journals.')
    a.receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True)


if __name__=='__main__':main()
