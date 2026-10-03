#!/usr/bin/env python3
"""Observe exact original joint-grid workers without changing or launching inference."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import time

from ancestral_chain_attempt import sha, write_json
from baliphy_sampler_resource_observation import describe, probe_native, read_cgroup, finish_roles, summarize, SUMMARY_FIELDS
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_reference_preflight import verify


PRODUCER_STATUS='complete_full_joint_sampler_resource_observation_pending_readback'
READER_STATUS='passed_full_joint_sampler_resource_observation_serialized_readback'
COMPLETED_STATUS='complete_verified_full_joint_sampler_resource_observation'


def events(path):
    with path.open() as handle:
        for i,line in enumerate(handle):
            value=json.loads(line);assert value['sequence']==i;yield value


def runtime_caps(reader=False):
    cpus,memory=(2,4) if reader else (1,2)
    group=next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    root=Path('/sys/fs/cgroup')/group.lstrip('/')
    limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':str(cpus*100000)+' 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
    assert all(os.environ.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    return limits


def replay(plan,jobs,root):
    """Stream every saved observation; bind actual attempt metadata at termination."""
    native=json.loads(Path(plan['sampler_plan']).read_text());record=json.loads(Path(plan['sampler_launch']).read_text())
    record['launch']=plan['sampler_launch'];assert fingerprint(record) is None
    terminal=json.loads((root/'sampler_terminal.json').read_text())
    actual=journal_terminal(record,expected_failure=terminal['status']=='verified_original_failure_preserved')
    for key in ['launch','pid','created','cmdline','status']:
        assert terminal[key]==actual[key]
    # Journal hashes can gain unrelated late manager messages; exact original
    # command/invocation and terminal resource/failure records stay mandatory.
    descriptors={};by_id={j['chain']['chain_id']:j for j in jobs}
    journal=root/'observations.jsonl';bindings=dict(plan['pins'])
    for snapshot in events(journal):
        assert snapshot['sampler_controller_pid']==record['pid'] and snapshot['sampler_controller_created']==record['created']
        assert snapshot['elapsed_since_previous_snapshot_seconds']>=0
        assert snapshot['cgroup']['limits']==record['actual_cgroup_limits']
        for descriptor in snapshot['new_attempts']:
            cid=descriptor['chain_id'];assert cid in by_id and cid not in descriptors
            assert descriptor==describe(by_id[cid],descriptor['identity_path']);descriptors[cid]=descriptor
        for sample in snapshot['native_observations']:
            cid=sample['chain_id'];assert cid in descriptors
            if sample['status']=='verified_live_native_resource_observation':
                job=by_id[cid];expected={'Max address space':job['memory_reservation_bytes'],
                    'Max cpu time':job['config']['timeout_seconds'],'Max file size':2*2**30}
                assert sample['actual_limits']==expected and sample['cgroup']==snapshot['cgroup']['path']
                assert set(sample['reported_memory_bytes'])=={'VmPeak','VmSize','VmHWM','VmRSS','VmSwap'}
    terminal_attempts=root/'attempts_at_termination.json'
    for descriptor in json.loads(terminal_attempts.read_text()):
        cid=descriptor['chain_id'];assert cid in by_id
        assert descriptor==describe(by_id[cid],descriptor['identity_path'])
        assert cid not in descriptors or descriptors[cid]==descriptor
        descriptors[cid]=descriptor
    receipts={}
    for cid,descriptor in descriptors.items():
        path=Path(descriptor['identity_path']).parent/'receipt.json'
        for field in ['identity','configuration','command']:
            bindings[descriptor[field+'_path']]=descriptor[field+'_sha256']
        if path.exists():receipts[cid]=dict(path=str(path),sha256=sha(path));bindings[str(path)]=sha(path)
    # A role with no live read remains an explicit missing observation, even
    # when its original native receipt proves that an attempt finished.
    def samples():
        for snapshot in events(journal):yield from snapshot['native_observations']
    rows=finish_roles(jobs,descriptors,samples(),receipts)
    summary=summarize(rows,events(journal),terminal['status'])
    source=Path(native['output'])/'receipt.json'
    if terminal['status']=='verified_original_terminal_success_with_bound_completed_artifacts':
        completed=json.loads(source.read_text());assert completed['full_chains']==1620
        assert completed['plan_sha256']==sha(plan['sampler_plan'])
        assert summary['roles_with_native_attempt']==1620 and len(receipts)==1620
        bindings[str(source)]=sha(source)
    for p in [journal,root/'sampler_terminal.json',terminal_attempts]:bindings[str(p)]=sha(p)
    verify(dict(pins=bindings));return rows,summary,bindings


def capture(plan,jobs,root):
    native=json.loads(Path(plan['sampler_plan']).read_text());verify(native)
    record=json.loads(Path(plan['sampler_launch']).read_text());record['launch']=plan['sampler_launch']
    assert record['actual_cgroup_limits']=={'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'}
    by_id={j['chain']['chain_id']:j for j in jobs};assert len(by_id)==1620
    journal=root/'observations.jsonl';assert not journal.exists(),'Partial telemetry requires review, never silent resumption'
    seen=set();sequence=0;previous=time.monotonic();last_report=previous
    with journal.open('x') as handle:
        while True:
            original=fingerprint(record)
            if original is None:
                try:terminal=journal_terminal(record)
                except AssertionError:terminal=journal_terminal(record,expected_failure=True)
                write_json(root/'sampler_terminal.json',terminal)
                # Fast attempts can start and finish between polling intervals.
                # Their receipts/identities still enter full role accounting;
                # never invent a live memory observation for them.
                descriptors=[describe(by_id[p.parent.parent.name],p) for p in
                    sorted((Path(native['output'])/'attempts').glob('*/attempt-*/process.json'))]
                assert len(descriptors)==len({d['chain_id'] for d in descriptors})
                write_json(root/'attempts_at_termination.json',descriptors);break
            assert shutil.disk_usage(root).free>=plan['resources']['minimum_free_disk_gib']*2**30
            try:cgroup=read_cgroup(record)
            except FileNotFoundError:
                # Original controller/cgroup may disappear between two reads;
                # repoll identity and actual terminal journals next iteration.
                time.sleep(plan['resources']['poll_seconds']);continue
            new=[];observed=[]
            for path in sorted((Path(native['output'])/'attempts').glob('*/attempt-*/process.json')):
                cid=path.parent.parent.name;assert cid in by_id
                if path not in seen:
                    assert not any(x['chain_id']==cid for x in new),'Multiple native attempts for one role'
                    descriptor=describe(by_id[cid],path);new.append(descriptor);seen.add(path)
                if not (path.parent/'receipt.json').exists():
                    descriptor=describe(by_id[cid],path)
                    observed.append(probe_native(by_id[cid],descriptor,cgroup['path']))
            # These sequential reads are not a simultaneous concurrency proof.
            now=time.monotonic();snapshot=dict(sequence=sequence,checked_utc=datetime.now(timezone.utc).isoformat(),
                sampler_controller_pid=record['pid'],sampler_controller_created=record['created'],
                elapsed_since_previous_snapshot_seconds=now-previous,cgroup=cgroup,new_attempts=new,native_observations=observed)
            handle.write(json.dumps(snapshot,allow_nan=False)+'\n');handle.flush()
            sequence+=1;previous=now
            if now-last_report>=60:
                print('sampler_resource_snapshots',sequence,'original_native_attempts_seen',len(seen),flush=True);last_report=now
            time.sleep(plan['resources']['poll_seconds'])


def run(path,reader=False):
    plan=json.loads(path.read_text());verify(plan);digest=sha(path);caps=runtime_caps(reader)
    native=json.loads(Path(plan['sampler_plan']).read_text());jobs=json.loads(Path(native['jobs']).read_text())
    root=Path(plan['output']);root.mkdir(exist_ok=True)
    lock=(root/('readback.lock' if reader else 'stage.lock')).open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    final=root/('readback.json' if reader else 'receipt.json');assert not final.exists(),'Completed telemetry cannot restart'
    marker=root/'stage_plan.json'
    if reader:assert json.loads(marker.read_text())==plan
    else:
        if marker.exists():assert json.loads(marker.read_text())==plan
        else:write_json(marker,plan)
        capture(plan,jobs,root)
    rows,summary,bindings=replay(plan,jobs,root);bindings[str(path)]=digest
    exported=root/'roles.json';result=dict(status=READER_STATUS if reader else PRODUCER_STATUS,
        plan_sha256=digest,**summary,source_hashes=bindings,actual_cgroup_limits=caps,
        scientific_eligibility=False,scope=plan['scope'])
    if reader:
        producer_path=root/'receipt.json';producer=json.loads(producer_path.read_text())
        assert producer['status']==PRODUCER_STATUS and producer['plan_sha256']==digest
        assert json.loads(exported.read_text())==rows
        assert all(producer[k]==summary[k] for k in SUMMARY_FIELDS)
        for name,h in producer['artifacts'].items():assert sha(root/name)==h
        result['producer_receipt_sha256']=sha(producer_path)
    else:
        write_json(exported,rows)
        result['artifacts']={p.name:sha(p) for p in [marker,exported,root/'observations.jsonl',root/'sampler_terminal.json',root/'attempts_at_termination.json']}
    verify(dict(pins=bindings));assert sha(path)==digest
    with final.open('x') as handle:handle.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--reader',action='store_true');args=parser.parse_args();run(args.plan,args.reader)
