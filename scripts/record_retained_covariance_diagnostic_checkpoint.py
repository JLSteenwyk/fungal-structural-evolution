#!/usr/bin/env python3
"""Observe the exact original live covariance diagnostic without restarting it."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    ep=Path('metadata/retained_timing_covariance_diagnostic_execution_20261004_v1.json')
    folder=ep.with_suffix('');config_path=folder/'configuration.json';child_path=folder/'process.json'
    config=json.loads(config_path.read_text());child=json.loads(child_path.read_text())
    pins=dict(config['source_hashes'])
    pp=Path('metadata/retained_timing_covariance_diagnostic_plan_20261004_v1.json')
    plan=json.loads(pp.read_text())
    for n,h in plan['pins'].items():bind(pins,n,h)
    for path in [pp,config_path,child_path,Path(__file__)]:bind(pins,path)
    verify(pins)
    wrapper=psutil.Process(config['wrapper']['pid'])
    assert wrapper.create_time()==config['wrapper']['created'] and wrapper.cmdline()==config['wrapper']['cmdline']
    native=psutil.Process(child['pid']);assert native.create_time()==child['created']
    assert native.cmdline()==child['command'][child['command'].index('--')+1:]
    assert native.status()!=psutil.STATUS_ZOMBIE
    group=next(l[3:] for l in Path(f'/proc/{native.pid}/cgroup').read_text().splitlines() if l.startswith('0::'))
    cg=Path('/sys/fs/cgroup')/group.lstrip('/')
    limits={k:(cg/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits==config['actual_cgroup_limits']
    unit='fungal-retained-timing-covariance-diagnostic-20261004-v1.service'
    rows=[json.loads(line) for line in subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True).splitlines()]
    inv=config['invocation_id'];rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in rows if r.get('_PID')==str(wrapper.pid) and r.get('_CMDLINE')==' '.join(wrapper.cmdline())]
    assert len(exact)==1 and json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=config['wrapper'],invocation_id=inv)
    assert sum(r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','') for r in rows)==1
    assert not any(r.get('CPU_USAGE_NSEC') for r in rows)
    root=Path(plan['output']);log=(folder/'stdout.log').read_text().splitlines()
    result=dict(status='verified_original_live_retained_covariance_diagnostic',checked_utc=datetime.now(timezone.utc).isoformat(),
        actual_original_tool_session_id=65481,actual_tool_wait_terminal=False,wrapper=config['wrapper'],
        invocation_id=inv,child=dict(pid=native.pid,created=native.create_time(),cmdline=native.cmdline(),
            status=native.status(),cpu_seconds=sum(native.cpu_times()[:2]),rss_bytes=native.memory_info().rss),
        actual_cgroup_limits=limits,memory_events={k:int(v) for k,v in (line.split() for line in (cg/'memory.events').read_text().splitlines())},
        last_progress_lines=log[-8:],completed_group_records=len(list((root/'groups').glob('*.json'))),
        receipt_present=Path('metadata/retained_timing_covariance_diagnostic_20261004_v1.json').exists(),
        source_hashes=pins,original_jobs_restarted=False,scientific_eligibility=False,gpu=False,
        scope='Fresh exact live original wrapper/native PID/create/command/invocation and actual limits. '
              'Snapshot of original source-loading/group progress only; no terminal wait, diagnostic '
              'completion, mismatch repair, original stage restart or biological acceptance.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','wrapper']},indent=2))


if __name__=='__main__':main()
