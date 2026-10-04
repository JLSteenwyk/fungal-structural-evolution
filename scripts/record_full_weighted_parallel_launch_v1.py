#!/usr/bin/env python3
"""Bind exact live original controller, invocation and actual resource limits."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--unit',required=True)
    p.add_argument('--wait-plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',type=int,required=True);p.add_argument('--memory-gib',type=int,required=True)
    p.add_argument('--original-tool-session',type=int,required=True)
    a=p.parse_args();assert not a.output.exists()
    pid=int(subprocess.check_output(['systemctl','--user','show',a.unit,'-p','MainPID','--value'],text=True))
    inv=subprocess.check_output(['systemctl','--user','show',a.unit,'-p','InvocationID','--value'],text=True).strip()
    assert pid>0 and len(inv)==32;process=psutil.Process(pid)
    command=[sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',str(a.wait_plan)]
    assert process.cmdline()==command and process.status()!=psutil.STATUS_ZOMBIE
    group=next(row[3:] for row in Path('/proc/'+str(pid)+'/cgroup').read_text().splitlines() if row.startswith('0::'))
    root=Path('/sys/fs/cgroup')/group.lstrip('/')
    limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':str(a.cpus*100000)+' 100000','memory.max':str(a.memory_gib*2**30),'memory.swap.max':'0'}
    assert all(process.environ()[key]=='1' for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    record=dict(unit=a.unit,pid=pid,created=process.create_time(),cmdline=command,plan=str(a.wait_plan),
        plan_sha256=sha(a.wait_plan),checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_cgroup_limits=limits,invocation_id=inv,original_tool_session_id=a.original_tool_session,
        source_hashes={str(Path(__file__)):sha(__file__)},scientific_eligibility=False)
    with a.output.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
