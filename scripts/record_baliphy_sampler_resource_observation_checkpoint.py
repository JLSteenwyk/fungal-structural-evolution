#!/usr/bin/env python3
"""Check exact original resource observer/controller handles and its latest snapshot."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal,live_record
from run_baliphy_reference_preflight import verify


def last_event(path):
    """Read a bounded suffix; ignore an unfinished line from a live writer."""
    with path.open('rb') as handle:
        handle.seek(0,2);size=handle.tell();handle.seek(max(0,size-2*2**20));tail=handle.read()
    lines=tail.split(b'\n')
    if size>2*2**20:lines=lines[1:]
    complete=[line for line in lines[:-1] if line]
    if not complete:return None
    return json.loads(complete[-1])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();path=Path('metadata/baliphy_sampler_resource_observation_plan_20261003.json')
    plan=json.loads(path.read_text());verify(plan);inventory=json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256']==sha(path);handles=[]
    launches=[*inventory['launches'],plan['sampler_launch']]
    for i,p in enumerate(launches):
        record=json.loads(Path(p).read_text());record['launch']=p
        assert sha(record['plan'])==record['plan_sha256'];proc=fingerprint(record)
        if proc is None:observed=journal_terminal(record)
        else:
            observed=live_record(record,proc)
            group=next(x[3:] for x in (Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines() if x.startswith('0::'))
            root=Path('/sys/fs/cgroup')/group.lstrip('/')
            limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            cpus,memory=[(1,2),(2,4),(2,4),(16,200)][i]
            assert limits==record['actual_cgroup_limits']=={'cpu.max':str(cpus*100000)+' 100000',
                'memory.max':str(memory*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits']=limits
        handles.append(observed)
    journal=Path(plan['output'])/'observations.jsonl';latest=last_event(journal) if journal.exists() else None
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),status='verified_original_sampler_resource_observer_runtime',
        plan_sha256=sha(path),frozen_pins_checked=len(plan['pins']),original_handles=handles,
        journal_bytes=journal.stat().st_size if journal.exists() else 0,
        latest_complete_snapshot=latest,posterior_qualified=False,new_native_jobs=0,gpu=False,
        scope='Exact original handles, live caps and latest complete observation only; no full event/native replay, source closure or precision/final-peak claim.')
    with args.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handles','latest_complete_snapshot']}))


if __name__=='__main__':main()
