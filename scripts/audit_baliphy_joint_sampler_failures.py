#!/usr/bin/env python3
"""Read-only completed-role failure census and immutable telemetry-prefix replay."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_qualification_v3 import inspect, SUCCESS
from baliphy_sampler_resource_observation import describe, empty_role, accumulate
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    plan_path=Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    plan=json.loads(plan_path.read_text());verify(plan['pins']);plan_digest=sha(plan_path)
    jobs=json.loads(Path(plan['jobs']).read_text());by={j['chain']['chain_id']:j for j in jobs}
    assert len(by)==len(jobs)==1620
    root=Path(plan['output']);paths=sorted((root/'chains').glob('*.json'))
    bindings={str(plan_path):plan_digest,str(Path(plan['jobs'])):sha(plan['jobs']),str(Path(__file__)):sha(__file__)}
    counts=Counter();failed=[];descriptors={};telemetry={}
    for path in paths:
        row=json.loads(path.read_text());counts[row['status']]+=1
        assert row['scientific_eligibility'] is row['posterior_qualified'] is False
        bindings[str(path)]=sha(path)
        if row['status']==SUCCESS:continue
        job=by[row['chain_id']];rp=Path(row['native_receipt']);receipt=json.loads(rp.read_text())
        assert row==inspect(job,rp,plan_digest,plan['mapping'],root/'frames'/row['chain_id'],allow_export_creation=False)
        assert sha(rp)==row['native_receipt_sha256'];bindings[str(rp)]=sha(rp)
        for name,h in receipt['artifacts'].items():
            fp=rp.parent/name;assert sha(fp)==h;bindings[str(fp)]=h
        descriptor=describe(job,rp.parent/'process.json');descriptors[row['chain_id']]=descriptor
        telemetry[row['chain_id']]=empty_role(job,descriptor)
        failed.append(dict(chain_id=row['chain_id'],family=row['family'],prior_label=row['prior_label'],
            chain_role=row['chain_role'],effective_input_group=row['effective_input_group'],
            proteins=job['chain']['proteins'],original_configuration_ids=row['original_configuration_ids'],
            status=row['status'],native_status=row['native_status'],exit_code=row['exit_code'],
            elapsed_worker_seconds=row['elapsed_worker_seconds'],native_receipt=str(rp),native_receipt_sha256=sha(rp),
            allocation_warning_lines=row['allocation_warning_lines'],bad_alloc=row['bad_alloc'],
            native_artifact_sizes={n:(rp.parent/n).stat().st_size for n in receipt['artifacts']},
            scientific_eligibility=False,posterior_qualified=False))
    observer_path=Path('metadata/baliphy_joint_sampler_resource_observation_v3_plan_20261003.json')
    observer=json.loads(observer_path.read_text());verify(observer['pins'])
    assert observer['sampler_plan']==str(plan_path)
    bindings[str(observer_path)]=sha(observer_path)
    launch=json.loads(Path(observer['sampler_launch']).read_text())
    original=Path(observer['output'])/'observations.jsonl';size=original.stat().st_size
    assert size<2**30,'Fresh larger bounded inventory required'
    snapshot=a.output/'observations-prefix.jsonl';snapshot_bytes=snapshots=observations=0
    first_events=last_events=None;maximum_peak=0;hasher=hashlib.sha256()
    with original.open('rb') as f,snapshot.open('xb') as g:
        while f.tell()<size:
            line=f.readline(size-f.tell())
            if not line.endswith(b'\n'):break
            event=json.loads(line)
            assert event['sequence']==snapshots
            assert event['sampler_controller_pid']==launch['pid'] and event['sampler_controller_created']==launch['created']
            assert event['cgroup']['limits']==launch['actual_cgroup_limits']
            if first_events is None:first_events=event['cgroup']['memory_events']
            last_events=event['cgroup']['memory_events']
            maximum_peak=max(maximum_peak,event['cgroup']['memory_bytes']['memory.peak'])
            for sample in event['native_observations']:
                if sample['chain_id'] not in telemetry:continue
                row=telemetry[sample['chain_id']];descriptor=descriptors[sample['chain_id']]
                assert sample['pid']==descriptor['pid'] and sample['created']==descriptor['created']
                if sample['status']=='verified_live_native_resource_observation':
                    assert sample['command']==by[sample['chain_id']]['config']['command'][5:]
                    assert sample['actual_limits']['Max address space']==by[sample['chain_id']]['memory_reservation_bytes']
                    assert sample['cgroup']==event['cgroup']['path']
                accumulate(telemetry,sample);observations+=1
            g.write(line);hasher.update(line);snapshot_bytes+=len(line);snapshots+=1
    assert sha(snapshot)==hasher.hexdigest()
    # A live journal can grow; the retained exact prefix remains reproducible.
    with original.open('rb') as f:
        h=hashlib.sha256();remaining=snapshot_bytes
        while remaining:
            block=f.read(min(1024*1024,remaining));assert block;h.update(block);remaining-=len(block)
        assert h.hexdigest()==hasher.hexdigest()
    bindings[str(snapshot)]=sha(snapshot);verify(bindings)
    for row in failed:row['resource_observations']=telemetry[row['chain_id']]
    result=dict(status='completed_read_only_joint_sampler_failure_and_resource_prefix_audit',
        checked_utc=datetime.now(timezone.utc).isoformat(),expected_roles=1620,completed_role_snapshot=len(paths),
        completed_role_status_counts=dict(counts),failed_or_invalid_roles=len(failed),failed_role_rows=failed,
        failed_exit_code_counts=dict(Counter(str(r['exit_code']) for r in failed)),
        failed_input_group_counts=dict(Counter(r['effective_input_group'] for r in failed)),
        original_observer_journal=str(original),original_journal_bytes_at_snapshot=size,
        immutable_telemetry_prefix=str(snapshot),immutable_prefix_bytes=snapshot_bytes,
        immutable_prefix_sha256=sha(snapshot),telemetry_snapshots=snapshots,failed_role_observations=observations,
        cgroup_initial_memory_events=first_events,cgroup_final_memory_events=last_events,
        maximum_sampled_cgroup_reported_peak_bytes=maximum_peak,source_hashes=bindings,
        existing_jobs_restarted=False,new_native_attempts=0,posterior_qualified=False,
        scope='Every completed role disposition at snapshot accounted; every failed/invalid row reconstructed from original native output and all its artifacts rehashed. Successful rows are counted,not freshly fully replayed here. Every complete telemetry line in an immutable live-journal prefix is replayed for failed original PID/create/command identities. Sampled Linux peaks can miss later peaks; a live prefix is not full observer closure or a memory guarantee. No source edits,acceptance,retry or posterior samples.')
    with a.receipt.open('x') as f:f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','failed_role_rows']},indent=2))


if __name__=='__main__':main()
