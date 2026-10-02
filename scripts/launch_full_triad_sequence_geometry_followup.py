#!/usr/bin/env python3
"""Queue frozen full preflight/geometry stages behind original native closure."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil
from run_ortholog_pair_guide_comparison import sha


def create(path,value):
    with Path(path).open('x') as f:f.write(json.dumps(value,indent=2)+'\n')


def launch(name,wait_path,command,dependencies,source_plan):
    unit='fungal-'+name+'-20261002-v1.service'
    launch_path='metadata/'+name.replace('-','_')+'_launch_20261002.json'
    assert not Path(launch_path).exists(),'Existing launch must remain immutable'
    pins={str(p):sha(p) for p in [source_plan,*dependencies,'scripts/run_after_verified_dependencies_v2.py']}
    create(wait_path,dict(dependencies=dependencies,command=command,pins=pins))
    cmd=['systemd-run','--user','--collect','--unit='+unit,'--working-directory='+str(Path.cwd()),
        '--property=CPUQuota=200%','--property=MemoryMax=32G','--property=MemorySwapMax=0',
        '--setenv=PYTHONUNBUFFERED=1','--setenv=OPENBLAS_NUM_THREADS=1','--setenv=OMP_NUM_THREADS=1','--setenv=MKL_NUM_THREADS=1',
        sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',str(wait_path)]
    subprocess.run(cmd,check=True)
    pid=0
    for _ in range(20):
        pid=int(subprocess.check_output(['systemctl','--user','show',unit,'-p','MainPID','--value'],text=True))
        if pid:break
        time.sleep(.1)
    assert pid>0
    proc=psutil.Process(pid);expected=[sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',str(wait_path)]
    assert proc.cmdline()==expected
    cg=subprocess.check_output(['systemctl','--user','show',unit,'-p','ControlGroup','--value'],text=True).strip()
    limits={name:(Path('/sys/fs/cgroup')/cg.lstrip('/')/name).read_text().strip() for name in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':'200000 100000','memory.max':str(32*2**30),'memory.swap.max':'0'}
    record=dict(unit=unit,pid=pid,created=proc.create_time(),cmdline=proc.cmdline(),plan=str(wait_path),plan_sha256=sha(wait_path),
        checked_utc=datetime.now(timezone.utc).isoformat(),actual_cgroup_limits=limits,systemd_launch_command=cmd)
    create(launch_path,record)
    print(json.dumps(dict(launch=launch_path,pid=pid,unit=unit)),flush=True)
    return launch_path


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workflow',type=Path,required=True);args=p.parse_args()
    workflow=json.loads(args.workflow.read_text())
    for path,digest in workflow['pins'].items():assert sha(path)==digest
    dependency=workflow['native_closure_launch']
    original=json.loads(Path(dependency).read_text());proc=psutil.Process(original['pid'])
    assert proc.create_time()==original['created'] and proc.cmdline()==original['cmdline'],'Require original live closure handle at queue creation'
    handles=[]
    for stage in workflow['stages']:
        name=stage['name'];source=stage['source_plan'];root=Path(json.loads(Path(source).read_text())['output'])
        producer=launch(name,'metadata/'+name.replace('-','_')+'_wait_plan_20261002.json',
            [sys.executable,stage['producer_script'],'--plan',source],[dependency],source)
        reader=launch(name+'-readback','metadata/'+name.replace('-','_')+'_readback_wait_plan_20261002.json',
            [sys.executable,stage['reader_script'],'--plan',source,'--output',str(root/'readback.json')],[producer],source)
        close_plan=stage['completion_plan'];spec=dict(source_plan=source,producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
            producer_status=stage['producer_status'],reader_status=stage['reader_status'],completed_status=stage['completed_status'],summary_fields=stage['summary_fields'],
            launches=[producer,reader],pins={str(path):sha(path) for path in [source,producer,reader,
                'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=stage['completion_output'],scope=stage['scope'])
        create(close_plan,spec)
        dependency=launch(name+'-closure','metadata/'+name.replace('-','_')+'_closure_wait_plan_20261002.json',
            [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',close_plan],[producer,reader],close_plan)
        handles.extend([producer,reader,dependency])
    create(workflow['launch_inventory'],dict(status='queued_full_triad_sequence_fit_preflight_and_geometry',
        workflow_sha256=sha(args.workflow),launches=handles,checked_utc=datetime.now(timezone.utc).isoformat(),
        scope=workflow['scope']))


if __name__=='__main__':main()
