#!/usr/bin/env python3
"""Run bounded retained covariance process-source software qualification with original execution evidence."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

import psutil

from audit_selected_taxon_identity_snapshot_v2 import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resources',type=Path,required=True)
    parser.add_argument('--execution',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--stage',choices=['software','producer','readback','summary'],required=True)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    a=parser.parse_args()
    assert not a.execution.exists() and not a.receipt.exists()
    plan=json.loads(a.resources.read_text())
    assert psutil.virtual_memory().available>=plan['minimum_available_ram_gib']*2**30
    assert psutil.disk_usage('.').free>=plan['minimum_free_disk_gib']*2**30
    cg=next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    cgroot=Path('/sys/fs/cgroup')/cg.lstrip('/')
    limits={key:(cgroot/key).read_text().strip() for key in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':'200000 100000','memory.max':str(16*2**30),'memory.swap.max':'0'}
    assert all(os.environ.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    raw=a.command[1:] if a.command[0]=='--' else a.command
    assert raw[0]==sys.executable and Path(raw[1]).is_file()
    pins={str(p):sha(p) for p in [Path(__file__),a.resources,Path(raw[1]),Path(sys.executable),Path('/usr/bin/prlimit')]}
    command=['/usr/bin/prlimit','--as='+str(plan['address_space_gib']*2**30),
             '--cpu='+str(plan['cpu_seconds_per_stage']),'--fsize='+str(plan['per_file_limit_mib']*2**20),'--',*raw]
    root=a.execution.with_suffix('');root.mkdir(parents=True,exist_ok=False)
    wrapper=psutil.Process()
    config=dict(stage=a.stage,command=command,working_directory=str(Path.cwd()),source_hashes=pins,
                invocation_id=os.environ.get('INVOCATION_ID'),wrapper=dict(pid=wrapper.pid,created=wrapper.create_time(),cmdline=wrapper.cmdline()),
                actual_cgroup_limits=limits,timeout_seconds=plan['wall_seconds_per_stage'],started_utc=datetime.now(timezone.utc).isoformat())
    with (root/'configuration.json').open('x') as f:f.write(json.dumps(config,indent=2)+'\n')
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'stdout.log').open('x') as out,(root/'stderr.log').open('x') as err:
        proc=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True)
        child=psutil.Process(proc.pid)
        identity=dict(pid=proc.pid,created=child.create_time(),command=command)
        with (root/'process.json').open('x') as f:f.write(json.dumps(identity,indent=2)+'\n')
        timed_out=False
        try:code=proc.wait(timeout=plan['wall_seconds_per_stage'])
        except subprocess.TimeoutExpired:
            timed_out=True
            try:os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            code=proc.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    for p,digest in pins.items():assert sha(p)==digest,p
    receipt_sha=sha(a.receipt) if a.receipt.exists() else None
    result=dict(**config,status='exited_zero_with_receipt' if code==0 and receipt_sha else 'failed_stage_retained',
                checked_utc=datetime.now(timezone.utc).isoformat(),child=identity,exit_code=code,timed_out=timed_out,
                wall_seconds=time.monotonic()-start,child_cpu_seconds=(after.ru_utime-before.ru_utime)+(after.ru_stime-before.ru_stime),
                child_peak_rss_bytes=after.ru_maxrss*1024,receipt=str(a.receipt),receipt_sha256=receipt_sha,
                artifacts={str(p):sha(p) for p in root.iterdir()},
                scope='One original CPU-only retained covariance process-source software qualification with actual cgroup caps, prlimit command, PID/create/exit and retained logs. Child RSS is not cgroup peak. Scope and counts are established by the software checker receipt. The checker covers600synthetic audits1200links and explicit synthetic source closure/23malformed states. It does not replay8680real certificates or launch real source analyses. Its own exit is not full-grid arithmetic closure. No original-output edits, failed-job restart, GPU work, new charges or biological acceptance.')
    with a.execution.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','scope','command','wrapper','child']}),flush=True)
    if code!=0 or not receipt_sha:raise SystemExit(code if code>0 else 1)


if __name__=='__main__':main()
