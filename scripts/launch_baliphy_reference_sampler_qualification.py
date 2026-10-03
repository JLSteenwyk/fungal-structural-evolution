#!/usr/bin/env python3
"""Queue all computational qualification roles behind original startup closure."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil

from ancestral_chain_attempt import sha
from baliphy_reference_sampler_qualification import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_reference_preflight import verify
from run_baliphy_reference_sampler_qualification import PRODUCER_STATUS, READER_STATUS, COMPLETED_STATUS


def launch(name, command, dependencies, source, cpus=2, memory=32):
    prefix='metadata/'+name.replace('-', '_'); unit='fungal-'+name+'-20261003-v1.service'
    wait=prefix+'_wait_plan_20261003.json'; record=prefix+'_launch_20261003.json'
    assert not Path(record).exists()
    create(wait, dict(dependencies=dependencies, command=command,
        pins={p:sha(p) for p in [source, *dependencies, 'scripts/run_after_verified_dependencies_v2.py']}))
    cmd=['systemd-run','--user','--collect','--unit='+unit,'--working-directory='+str(Path.cwd()),
        '--property=CPUQuota='+str(cpus*100)+'%','--property=MemoryMax='+str(memory)+'G',
        '--property=MemorySwapMax=0','--setenv=PYTHONUNBUFFERED=1','--setenv=OPENBLAS_NUM_THREADS=1',
        '--setenv=OMP_NUM_THREADS=1','--setenv=MKL_NUM_THREADS=1',
        sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',wait]
    subprocess.run(cmd,check=True); pid=0
    for _ in range(20):
        pid=int(subprocess.check_output(['systemctl','--user','show',unit,'-p','MainPID','--value'],text=True))
        if pid:break
        time.sleep(.1)
    assert pid>0; process=psutil.Process(pid)
    expected=[sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',wait]
    assert process.cmdline()==expected
    group=subprocess.check_output(['systemctl','--user','show',unit,'-p','ControlGroup','--value'],text=True).strip()
    root=Path('/sys/fs/cgroup')/group.lstrip('/')
    limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':str(cpus*100000)+' 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
    create(record,dict(unit=unit,pid=pid,created=process.create_time(),cmdline=expected,plan=wait,
        plan_sha256=sha(wait),checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits,systemd_launch_command=cmd))
    print(json.dumps(dict(launch=record,pid=pid,unit=unit)),flush=True);return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());verify(plan)
    assert psutil.virtual_memory().available>=200*2**30 and shutil.disk_usage('.').free>=128*2**30
    assert not Path(plan['output']).exists()
    for path in plan['dependencies']:
        record=json.loads(Path(path).read_text());record['launch']=path
        if fingerprint(record) is None:journal_terminal(record)
    producer=launch('baliphy-reference-sampler-qualification',
        [sys.executable,'scripts/run_baliphy_reference_sampler_qualification.py','--plan',str(args.plan)],
        plan['dependencies'],str(args.plan),cpus=16,memory=200)
    reader=launch('baliphy-reference-sampler-qualification-readback',
        [sys.executable,'scripts/run_baliphy_reference_sampler_qualification.py','--plan',str(args.plan),'--reader'],
        [producer],str(args.plan))
    root=Path(plan['output']);completion=dict(source_plan=str(args.plan),producer_receipt=str(root/'receipt.json'),
        independent_readback=str(root/'readback.json'),producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,
        completed_status=COMPLETED_STATUS,summary_fields=SUMMARY_FIELDS+['reservation_audit'],launches=[producer,reader],
        pins={p:sha(p) for p in [str(args.plan),producer,reader,'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=plan['scope'])
    create(plan['completion_plan'],completion)
    closer=launch('baliphy-reference-sampler-qualification-closure',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['completion_plan']],
        [producer,reader],plan['completion_plan'])
    create(plan['launch_inventory'],dict(status='queued_full_reference_short_sampler_qualification',
        source_plan=str(args.plan),source_plan_sha256=sha(args.plan),launches=[producer,reader,closer],
        expected_native_roles=1620,iterations_per_role=20,posterior_qualified=False,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
