#!/usr/bin/env python3
"""Run new stack-limited attempts and retain complete 1620-role provenance."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import threading
import time

import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_sampler_qualification_v3 import inspect
from baliphy_sampler_resource_observation import describe, probe_native, read_cgroup, finish_roles
from baliphy_stack_followup import prerequisites, merge
from reference_measurement_union_sources import verify
from reference_sampler_memory_budget import MemoryBudget
from run_baliphy_joint_sampler_qualification_v3 import execute_job, verify_ledger, runtime_caps


PRODUCER_STATUS = 'complete_full_role_stack_followup_dispositions_pending_readback_v1'
READER_STATUS = 'passed_full_role_stack_followup_serialized_readback_v1'
COMPLETED_STATUS = 'complete_verified_full_role_stack_followup_v1'


def stack_limit(text):
    matches = re.findall(r'^Max stack size\s+(\d+)\s+(\d+)\s+bytes\s*$',text,re.M)
    assert len(matches) == 1 and tuple(map(int,matches[0])) == (64*2**20,64*2**20)
    return 64*2**20


def stack_probe(job,descriptor,group):
    sample = probe_native(job,descriptor,group)
    if sample['status'] != 'verified_live_native_resource_observation': return sample
    try:
        process = psutil.Process(descriptor['pid'])
        limits = (Path('/proc')/str(process.pid)/'limits').read_text()
        native = job['config']['command'][job['config']['command'].index('--')+1:]
        if process.create_time() != descriptor['created'] or process.status() == psutil.STATUS_ZOMBIE or process.cmdline() != native:
            return dict(chain_id=descriptor['chain_id'],pid=descriptor['pid'],created=descriptor['created'],
                        identity_sha256=descriptor['identity_sha256'],status='original_process_unavailable_after_observation')
        sample['native_stack_bytes'] = stack_limit(limits)
        return sample
    except (psutil.NoSuchProcess,psutil.ZombieProcess,ProcessLookupError,FileNotFoundError):
        return dict(chain_id=descriptor['chain_id'],pid=descriptor['pid'],created=descriptor['created'],
                    identity_sha256=descriptor['identity_sha256'],status='process_unavailable_during_observation')


def capture(root,jobs,controller,stop):
    by = {j['chain']['chain_id']:j for j in jobs}; sequence = 0; previous = time.monotonic()
    with (root/'observations.jsonl').open('x') as handle:
        while True:
            group = read_cgroup(controller); samples = []
            for path in sorted((root/'attempts').glob('*/attempt-0001/process.json')):
                cid = path.parent.parent.name; assert cid in by
                if not (path.parent/'receipt.json').exists():
                    samples.append(stack_probe(by[cid],describe(by[cid],path),group['path']))
            now = time.monotonic()
            event = dict(sequence=sequence,checked_utc=datetime.now(timezone.utc).isoformat(),
                controller_pid=controller['pid'],controller_created=controller['created'],
                elapsed_since_previous_snapshot_seconds=now-previous,cgroup=group,native_observations=samples)
            handle.write(json.dumps(event,allow_nan=False)+'\n'); handle.flush()
            sequence += 1; previous = now
            if stop.is_set(): break
            stop.wait(1)


def resource_replay(root,jobs,expected_caps):
    controller = json.loads((root/'controller.json').read_text())
    assert controller['actual_cgroup_limits'] == expected_caps
    by = {j['chain']['chain_id']:j for j in jobs}; descriptors = {}; receipts = {}; bindings = {}
    for job in jobs:
        cid = job['chain']['chain_id']; folder = root/'attempts'/cid
        assert {p.name for p in folder.glob('attempt-[0-9]*')} == {'attempt-0001'}
        ip = folder/'attempt-0001/process.json'; d = describe(job,ip); descriptors[cid] = d
        rp = ip.parent/'receipt.json'; receipts[cid] = dict(path=str(rp),sha256=sha(rp))
        for name in ['identity','configuration','command']:
            bindings[d[name+'_path']] = d[name+'_sha256']
        bindings[str(rp)] = sha(rp)
    count = 0; peak = maximum_gap = 0; samples = []
    for i,line in enumerate((root/'observations.jsonl').open()):
        event = json.loads(line); assert event['sequence'] == i
        assert event['controller_pid'] == controller['pid'] and event['controller_created'] == controller['created']
        group = event['cgroup']; assert group['limits'] == expected_caps
        assert group['memory_bytes']['memory.swap.current'] == 0
        assert all(group['memory_events'][k] == 0 for k in ['max','oom','oom_kill','oom_group_kill'])
        peak = max(peak,group['memory_bytes']['memory.peak'])
        assert peak <= 200*2**30
        maximum_gap = max(maximum_gap,event['elapsed_since_previous_snapshot_seconds']); count += 1
        for sample in event['native_observations']:
            cid = sample['chain_id']; d = descriptors[cid]; job = by[cid]
            assert sample['pid'] == d['pid'] and sample['created'] == d['created'] and sample['identity_sha256'] == d['identity_sha256']
            if sample['status'] == 'verified_live_native_resource_observation':
                assert sample['command'] == job['config']['command'][job['config']['command'].index('--')+1:]
                assert sample['actual_limits'] == {'Max address space':48*2**30,
                    'Max cpu time':job['config']['timeout_seconds'],'Max file size':2*2**30}
                assert sample['native_stack_bytes'] == 64*2**20 and sample['cgroup'] == group['path']
            samples.append(sample)
    assert count > 0
    roles = finish_roles(jobs,descriptors,samples,receipts)
    assert len(roles) == 24 and all(r['native_receipt'] is not None for r in roles)
    bindings[str(root/'controller.json')] = sha(root/'controller.json')
    bindings[str(root/'observations.jsonl')] = sha(root/'observations.jsonl')
    return roles,dict(expected_followup_roles=24,observations=count,
        roles_with_live_stack_limit_observation=sum(r['live_observations'] > 0 for r in roles),
        roles_without_live_observation=sum(r['live_observations'] == 0 for r in roles),
        observed_cgroup_reported_peak_bytes=peak,maximum_observed_poll_gap_seconds=maximum_gap,
        no_observed_cgroup_oom_or_limit_event=True,stack_limit_bytes=64*2**20,
        scope='Identity-bound sampled procfs and cgroup resources; missed/transitioning observations retained. Not precise final native peaks or longer-chain resource bounds.'),bindings


def run(path,reader=False):
    plan = json.loads(path.read_text()); verify(plan['pins']); digest = sha(path)
    caps = runtime_caps(plan,reader); originals,original_rows,jobs,bindings = prerequisites(plan)
    bindings.update(plan['pins']); bindings[str(path)] = digest
    root = Path(plan['output']); root.mkdir(exist_ok=True)
    lock = (root/('readback.lock' if reader else 'stage.lock')).open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    final = root/('readback.json' if reader else 'receipt.json'); assert not final.exists(),'Completed follow-up cannot restart'
    source = json.loads(Path(plan['source_plan']).read_text()); started = time.monotonic()
    if reader:
        pp = root/'receipt.json'; producer = json.loads(pp.read_text())
        assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == digest
        verify(producer['source_hashes'])
        assert json.loads((root/'stage_plan.json').read_text()) == plan
        rows = []
        for job in jobs:
            cid = job['chain']['chain_id']; cp = root/'chains'/(cid+'.json'); saved = json.loads(cp.read_text())
            row = inspect(job,Path(saved['native_receipt']),digest,plan['mapping'],root/'frames'/cid,allow_export_creation=False)
            assert row == saved; rows.append(row)
        for name,h in producer['artifacts'].items(): assert sha(root/name) == h
    else:
        assert psutil.virtual_memory().available >= 200*2**30
        assert not (root/'stage_plan.json').exists(), 'Partial producer requires review; never silently resume'
        write_json(root/'stage_plan.json',plan); (root/'chains').mkdir(exist_ok=False)
        proc = psutil.Process()
        controller = dict(pid=proc.pid,created=proc.create_time(),cmdline=proc.cmdline(),
            invocation_id=os.environ['INVOCATION_ID'],actual_cgroup_limits=caps)
        write_json(root/'controller.json',controller)
        stop = threading.Event(); rows = []
        with (root/'memory_reservations.jsonl').open('x') as handle:
            def persist(event):handle.write(json.dumps(event)+'\n');handle.flush()
            budget = MemoryBudget(192*2**30,persist)
            with ThreadPoolExecutor(max_workers=1) as observer, ThreadPoolExecutor(max_workers=4) as pool:
                observation = observer.submit(capture,root,jobs,controller,stop)
                def execute(job):
                    assert not (root/'attempts'/job['chain']['chain_id']).exists(),'Follow-up gets exactly one new attempt'
                    if observation.done(): observation.result(); raise RuntimeError('Resource observer stopped early')
                    return execute_job(job,root,digest,plan['mapping'],budget,plan['resources']['minimum_free_disk_gib'])
                try:
                    for future in as_completed([pool.submit(execute,j) for j in jobs]):
                        rows.append(future.result());print('stack_followup',len(rows),'/24',rows[-1]['status'],flush=True)
                        if observation.done(): observation.result(); raise RuntimeError('Resource observer stopped early')
                finally:
                    stop.set(); observation.result()
            assert budget.aborted is None
    rows.sort(key=lambda r:r['chain_id']); assert len(rows) == 24
    assert {p.name for p in (root/'chains').glob('*.json')} == {r['chain_id']+'.json' for r in rows}
    ledger,summary = merge(originals,original_rows,jobs,rows,source['output'],root)
    memory_events = [json.loads(line) for line in (root/'memory_reservations.jsonl').read_text().splitlines()]
    summary['reservation_audit'] = verify_ledger(memory_events,jobs,192*2**30,4)
    roles,summary['resource_audit'],rb = resource_replay(root,jobs,{'cpu.max':'400000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'})
    bindings.update(rb)
    if reader:
        assert json.loads((root/'dispositions.json').read_text()) == rows
        assert json.loads((root/'full_role_accounting.json').read_text()) == ledger
        assert json.loads((root/'resources.json').read_text()) == roles
        from baliphy_stack_followup import SUMMARY_FIELDS
        assert all(producer[k] == summary[k] for k in SUMMARY_FIELDS)
    else:
        for name,value in [('dispositions.json',rows),('full_role_accounting.json',ledger),('resources.json',roles)]:
            write_json(root/name,value)
    for row in rows:
        rp = Path(row['native_receipt']); bindings[str(rp)] = sha(rp)
        for name,h in json.loads(rp.read_text())['artifacts'].items(): bindings[str(rp.parent/name)] = h
        for frame in row['joint_frames']: bindings[frame['projection_array']] = frame['projection_array_sha256']
    expected = {Path(f['projection_array']) for r in rows for f in r['joint_frames']}
    assert set((root/'frames').glob('*/*')) == expected
    verify(bindings)
    result = dict(status=READER_STATUS if reader else PRODUCER_STATUS,plan_sha256=digest,**summary,
        actual_cgroup_limits=caps,source_hashes=bindings,elapsed_seconds=time.monotonic()-started,
        scientific_eligibility=False,scope=plan['scope'])
    if reader:result['producer_receipt_sha256'] = sha(root/'receipt.json')
    else:result['artifacts'] = {str(p.relative_to(root)):sha(p) for p in [root/'stage_plan.json',root/'dispositions.json',
        root/'full_role_accounting.json',root/'resources.json',root/'observations.jsonl',root/'controller.json',
        root/'memory_reservations.jsonl',*sorted((root/'chains').glob('*.json'))]}
    with final.open('x') as handle:handle.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary),flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--reader',action='store_true');a = p.parse_args();run(a.plan,a.reader)
