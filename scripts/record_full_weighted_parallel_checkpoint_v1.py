#!/usr/bin/env python3
"""Observe exact original parallel controller, native parent and worker limits."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def identity(process):
    return dict(pid=process.pid,created=process.create_time(),cmdline=process.cmdline(),
        parent_pid=process.ppid(),status=process.status(),cpu_seconds=sum(process.cpu_times()[:2]),
        rss_bytes=process.memory_info().rss,vm_bytes=process.memory_info().vms,
        address_space_limits=list(process.rlimit(resource.RLIMIT_AS)))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    pp=Path('metadata/full_weighted_parallel_qualification_plan_20261004_v1.json')
    ip=Path('metadata/full_weighted_parallel_qualification_launches_20261004_v1.json')
    plan=json.loads(pp.read_text());inventory=json.loads(ip.read_text());verify(plan['pins'])
    assert inventory['source_plan_sha256']==sha(pp)
    handles=[]
    for index,path in enumerate(inventory['launches']):
        limits={'cpu.max':str((16 if index<2 else 2)*100000)+' 100000',
            'memory.max':str((200 if index<2 else 32)*2**30),'memory.swap.max':'0'}
        handles.append(observe(path,limits))
    launch=json.loads(Path(inventory['launches'][0]).read_text())
    controller=psutil.Process(launch['pid'])
    assert controller.create_time()==launch['created'] and controller.cmdline()==launch['cmdline']
    command=[sys.executable,'scripts/full_weighted_covariance_qualification_parallel_v1.py','--plan',str(pp)]
    parents=[p for p in controller.children() if p.cmdline()==command]
    assert len(parents)==1;parent=parents[0]
    assert parent.rlimit(resource.RLIMIT_AS)==(8*2**30,24*2**30)
    workers=parent.children();assert 0<len(workers)<=16
    for worker in workers:
        assert worker.cmdline()==command and worker.rlimit(resource.RLIMIT_AS)==(12*2**30,12*2**30)
        assert all(worker.environ()[k]=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    group=next(row[3:] for row in Path('/proc/'+str(parent.pid)+'/cgroup').read_text().splitlines() if row.startswith('0::'))
    cg=Path('/sys/fs/cgroup')/group.lstrip('/')
    memory={key:(cg/key).read_text().strip() for key in ['memory.current','memory.peak','memory.swap.current','memory.events']}
    root=Path(plan['output']);completed=[];partial=[];pins=dict(plan['pins'])
    for path in sorted((root/'checkpoints').glob('*.producer.json')):
        try:row=json.loads(path.read_text())
        except json.JSONDecodeError:partial.append(str(path));continue
        assert row['cached_numeric_inputs_preserved'] and row['worker_address_space_limit_bytes']==12*2**30
        assert row['audits']==1200 and not row['scientific_eligibility']
        completed.append(dict(cohort_index=row['cohort_index'],cohort_id=row['cohort_id'],
            audits=row['audits'],links=row['links'],worker=row['worker']))
        bind(pins,path)
    assert len({row['cohort_index'] for row in completed})==len(completed)
    for path in [Path(__file__),pp,ip,*map(Path,inventory['launches'])]:bind(pins,path)
    verify(pins)
    result=dict(status='verified_live_original_full_parallel_weighted_numerical_calculation',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_plan_sha256=sha(pp),
        original_handles=handles,native_parent=identity(parent),native_workers=[identity(p) for p in workers],
        worker_count=len(workers),actual_parent_soft_address_space_gib=8,actual_worker_address_space_gib=12,
        cgroup_memory=memory,completed_cohort_checkpoints=len(completed),completed_cohorts=completed,
        partially_serialized_checkpoint_files=partial,expected=plan['expected'],
        failure_files=[str(path) for path in (root/'failures').rglob('*') if path.is_file()],
        producer_receipt_present=(root/'receipt.json').exists(),independent_readback_present=(root/'readback.json').exists(),
        completion_present=Path(plan['completion']).exists(),source_hashes=pins,
        working_model_fits_computed=0,scientific_eligibility=False,gpu=False,
        scope='Exact live original controller/native parent/fork worker identities and actual '
              'limits are checked. Cohort checkpoints indicate completed producer arithmetic '
              'and preserved worker caches, not independent full-grid acceptance. Source hashes '
              'and available checkpoint hashes are reverified. Reader/closure, full timing and '
              'biological fitting remain required. RSS sums would double-count shared fork pages.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','native_workers','original_handles','completed_cohorts']},indent=2))


if __name__=='__main__':main()
