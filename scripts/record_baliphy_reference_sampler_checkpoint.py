#!/usr/bin/env python3
"""Observe the original qualification queue, native limits and retained outcomes."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import psutil

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_baliphy_reference_preflight import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();path=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    plan=json.loads(path.read_text());verify(plan)
    inventory=json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256']==sha(path)
    handles=[]
    for i,launch in enumerate(inventory['launches']):
        record=json.loads(Path(launch).read_text());record['launch']=launch
        assert sha(record['plan'])==record['plan_sha256']
        process=fingerprint(record)
        if process is None:observed=journal_terminal(record)
        else:
            observed=live_record(record,process)
            group=next(x[3:] for x in (Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines() if x.startswith('0::'))
            root=Path('/sys/fs/cgroup')/group.lstrip('/')
            limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            assert limits==record['actual_cgroup_limits']=={'cpu.max':('1600000' if i==0 else '200000')+' 100000',
                'memory.max':str((200 if i==0 else 32)*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits']=limits
            observed['current_cgroup_memory_bytes']={k:(root/k).read_text().strip() for k in ['memory.current','memory.peak']}
        handles.append(observed)
    jobs={j['chain']['chain_id']:j for j in json.loads(Path(plan['jobs']).read_text())}
    root=Path(plan['output']);checkpoints=[json.loads(p.read_text()) for p in sorted((root/'chains').glob('*.json'))]
    assert all(r['plan_sha256']==sha(path) and r['posterior_qualified'] is False for r in checkpoints)
    native=[];vanished=[]
    for identity in sorted((root/'attempts').glob('*/attempt-*/process.json')):
        if (identity.parent/'receipt.json').exists():continue
        saved=json.loads(identity.read_text());cid=identity.parent.parent.name;job=jobs[cid]
        try:
            proc=psutil.Process(saved['pid'])
            if proc.create_time()!=saved['created'] or proc.status()==psutil.STATUS_ZOMBIE:continue
            command=job['config']['command'];assert saved['command']==command
            assert proc.cmdline() in [command,command[command.index('--')+1:]]
            text=(Path('/proc')/str(proc.pid)/'limits').read_text()
            limits={}
            for name,value in [('Max address space',job['memory_reservation_bytes']),
                               ('Max cpu time',job['config']['timeout_seconds']),('Max file size',2*2**30)]:
                match=re.search(r'^'+re.escape(name)+r'\s+(\d+)\s+(\d+)\s+',text,re.M)
                assert match and int(match[1])==int(match[2])==value
                limits[name]=value
            native.append(dict(chain_id=cid,pid=proc.pid,created=saved['created'],command=proc.cmdline(),actual_limits=limits))
        except (psutil.NoSuchProcess,ProcessLookupError,FileNotFoundError):vanished.append(saved['pid'])
    assert len(native)<=16 and sum(x['actual_limits']['Max address space'] for x in native)<=192*2**30
    completion=None
    if Path(plan['completion']).exists():
        completion=json.loads(Path(plan['completion']).read_text())
        assert sha(completion['full_hash_archive'])==completion['full_hash_archive_sha256']
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),status='verified_original_reference_sampler_queue_and_native_limits',
        source_plan=str(path),source_plan_sha256=sha(path),frozen_pins_checked=len(plan['pins']),original_handles=handles,
        retained_checkpoints=len(checkpoints),status_counts=dict(Counter(r['status'] for r in checkpoints)),
        observed_native_workers=native,ephemeral_native_processes_disappeared_during_observation=vanished,
        accounting_closure=completion,posterior_qualified=False,gpu=False,
        scope='Original queue identities/caps and actual live native process limits checked; checkpoint counts are unclosed progress. Not a full native/output replay or posterior qualification.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handles','observed_native_workers','accounting_closure']}))


if __name__=='__main__':main()
