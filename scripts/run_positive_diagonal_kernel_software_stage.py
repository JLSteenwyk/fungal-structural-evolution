#!/usr/bin/env python3
"""Bound complete declared positive-diagonal kernel software qualification."""
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

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import runtime_caps


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--resources', type=Path, required=True)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args(); assert not a.execution.exists() and not a.receipt.exists()
    plan = json.loads(a.resources.read_text()); caps = runtime_caps()
    assert psutil.virtual_memory().available >= plan['minimum_available_ram_gib'] * 2**30
    assert psutil.disk_usage('.').free >= plan['minimum_free_disk_gib'] * 2**30
    raw = a.command[1:] if a.command[0] == '--' else a.command
    assert raw[0] == sys.executable and Path(raw[1]).is_file()
    pins = {str(q):sha(q) for q in [Path(__file__),a.resources,Path(raw[1]),Path(sys.executable),Path('/usr/bin/prlimit')]}
    command = ['/usr/bin/prlimit', '--as=' + str(plan['address_space_gib'] * 2**30),
        '--cpu=' + str(plan['cpu_seconds_per_stage']), '--fsize=' + str(plan['per_file_limit_mib'] * 2**20), '--', *raw]
    root = a.execution.with_suffix(''); root.mkdir(parents=True, exist_ok=False)
    wrapper = psutil.Process()
    config = dict(command=command, working_directory=str(Path.cwd()), source_hashes=pins,
        invocation_id=os.environ['INVOCATION_ID'], wrapper=dict(pid=wrapper.pid,created=wrapper.create_time(),cmdline=wrapper.cmdline()),
        actual_cgroup_limits=caps, timeout_seconds=plan['wall_seconds_per_stage'],started_utc=datetime.now(timezone.utc).isoformat())
    with (root/'configuration.json').open('x') as f:json.dump(config,f,indent=2);f.write('\n')
    print(json.dumps(dict(original_wrapper=config['wrapper'],invocation_id=config['invocation_id'])),flush=True)
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (root/'stdout.log').open('x') as out,(root/'stderr.log').open('x') as err:
        proc=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True)
        native=psutil.Process(proc.pid); identity=dict(pid=proc.pid,created=native.create_time(),command=command)
        with (root/'process.json').open('x') as f:json.dump(identity,f,indent=2);f.write('\n')
        timed_out=False
        try:code=proc.wait(timeout=plan['wall_seconds_per_stage'])
        except subprocess.TimeoutExpired:
            timed_out=True
            try:os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            code=proc.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    for q,h in pins.items():assert sha(q)==h
    rh=sha(a.receipt) if a.receipt.exists() else None
    result=dict(**config,status='exited_zero_with_receipt' if code==0 and rh else 'failed_stage_retained',
        checked_utc=datetime.now(timezone.utc).isoformat(),child=identity,exit_code=code,timed_out=timed_out,
        wall_seconds=time.monotonic()-start,child_cpu_seconds=after.ru_utime-before.ru_utime+after.ru_stime-before.ru_stime,
        child_peak_rss_bytes=after.ru_maxrss*1024,receipt=str(a.receipt),receipt_sha256=rh,
        artifacts={str(q):sha(q) for q in root.iterdir()},
        scope='Original bounded complete declared positive-diagonal latent/raw/REML kernel checker, exact wrapper/native identities, exit, caps and logs. Native RSS is not a cgroup peak. No real full-source weighted qualification, fitting, restart, posterior acceptance or GPU use.')
    with a.execution.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}),flush=True)
    if code!=0 or rh is None:raise SystemExit(code if code>0 else 1)


if __name__=='__main__':main()
