#!/usr/bin/env python3
"""Observe original failed-role comparison identities, native caps and the active lease ledger."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    pp=Path('metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json')
    cp=Path('metadata/baliphy_joint_fasta_v7_failure_grid_execution_20261004_v1/configuration.json')
    tool_path=Path('metadata/baliphy_joint_fasta_v7_failure_grid_original_tool_initial_20261004_v1.json')
    plan,configuration,tool=[json.loads(p.read_text()) for p in [pp,cp,tool_path]]
    assert tool['session_id']==72605 and configuration['invocation_id'] in tool['output']
    wrapper=configuration['wrapper'];process=psutil.Process(wrapper['pid'])
    assert process.create_time()==wrapper['created'] and process.cmdline()==wrapper['cmdline']
    group=next(line[3:] for line in Path('/proc',str(process.pid),'cgroup').read_text().splitlines() if line.startswith('0::'))
    cgroup=Path('/sys/fs/cgroup')/group.lstrip('/')
    limits={k:(cgroup/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits==configuration['actual_cgroup_limits']=={'cpu.max':'400000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'}
    jobs=json.loads(Path(plan['jobs']).read_text());expected={tuple(j['config']['command'][j['config']['command'].index('--')+1:]):j for j in jobs}
    children=[];parent=[]
    for child in process.children(recursive=True):
        try:
            created=child.create_time();command=child.cmdline()
            if 'scripts/run_baliphy_joint_fasta_v7_failure_grid_v1.py' in command:
                assert child.rlimit(psutil.RLIMIT_AS)==(8*2**30,192*2**30)
                parent.append(dict(pid=child.pid,created=created,cmdline=command,status=child.status(),
                                   address_space_limits=child.rlimit(psutil.RLIMIT_AS)))
            elif tuple(command) in expected:
                job=expected[tuple(command)];caps=job['config']['command']
                cpu=int(next(v.split('=')[1] for v in caps if v.startswith('--cpu=')))
                assert child.rlimit(psutil.RLIMIT_AS)==(48*2**30,)*2
                assert child.rlimit(psutil.RLIMIT_CPU)==(cpu,)*2
                assert child.rlimit(psutil.RLIMIT_FSIZE)==(2*2**30,)*2
                assert child.rlimit(psutil.RLIMIT_STACK)==(8*2**20,-1)
                assert all(child.environ().get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
                children.append(dict(chain_id=job['chain']['chain_id'],pid=child.pid,created=created,cmdline=command,
                    status=child.status(),cpu_seconds=child.cpu_times().user+child.cpu_times().system,
                    rss_bytes=child.memory_info().rss,address_space_limits=child.rlimit(psutil.RLIMIT_AS),
                    cpu_limits=child.rlimit(psutil.RLIMIT_CPU),file_limits=child.rlimit(psutil.RLIMIT_FSIZE),
                    stack_limits=child.rlimit(psutil.RLIMIT_STACK)))
            assert child.create_time()==created
        except (psutil.NoSuchProcess,FileNotFoundError):continue
    assert len(parent)==1 and len(children)<=4
    root=Path(plan['output']);ledger=root/'memory_reservations.jsonl'
    events=[json.loads(line) for line in ledger.read_text().splitlines()]
    active={};seen=set();released=set();total=0
    for i,event in enumerate(events):
        assert event['sequence']==i;identifier=event['identifier'];amount=event['amount']
        assert amount==48*2**30
        if event['action']=='acquire':
            assert identifier not in seen;seen.add(identifier);active[identifier]=amount;total+=amount
        else:
            assert event['action']=='release' and identifier not in released
            assert active.pop(identifier)==amount;released.add(identifier);total-=amount
        assert event['reserved_total']==total and event['active_roles']==len(active)
        assert 0<=total<=192*2**30 and len(active)<=4
    assert all(child['chain_id'] in active for child in children)
    result=dict(status='verified_original_live_all24_v7_failure_comparison_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=72605,
        invocation_id=configuration['invocation_id'],original_wrapper=wrapper,controller=parent[0],
        actual_native_processes=children,current_cgroup_limits=limits,
        cgroup_memory={k:int((cgroup/k).read_text()) for k in ['memory.current','memory.peak','memory.swap.current']},
        active_lease_roles=len(active),current_reserved_bytes=total,observed_lease_events=len(events),
        completed_role_checkpoints=len(list((root/'chains').glob('*.json'))),expected_roles=24,
        producer_receipt_present=(root/'receipt.json').exists(),reader_receipt_present=(root/'readback.json').exists(),
        source_hashes={str(q):sha(q) for q in [Path(__file__),pp,cp,tool_path,Path(plan['jobs'])]},
        all24_output_integrity_qualified=False,scientific_eligibility=False,posterior_qualified=False,gpu=False,
        scope='Exact original live wrapper/controller/native identities,actual caps/BLAS settings '
              'and current lease-prefix accounting. Full native/source/output/independent/journal '
              'closure is not inferred from this observation. AS leases are not measured memory usage.')
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_wrapper','controller','actual_native_processes','source_hashes']},indent=2))


if __name__=='__main__':main()
