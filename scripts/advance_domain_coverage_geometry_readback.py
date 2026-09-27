#!/usr/bin/env python3
"""Wait for the exact joined-cohort producer and verify its complete output."""
import json,time,subprocess,sys,hashlib
from pathlib import Path
import psutil
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lp=Path('metadata/domain_coverage_geometry_launch_20260927.json');lh=sha(lp);dep=json.loads(lp.read_text());pp=Path(dep['cmdline'][-1]);ph=sha(pp);assert ph==dep['plan_sha256']
plan=json.loads(pp.read_text())
while True:
    try:
        p=psutil.Process(dep['pid'])
        if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
        assert p.cmdline()==dep['cmdline']
    except psutil.NoSuchProcess:break
    time.sleep(30)
assert sha(lp)==lh and sha(pp)==ph
state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
for p,h in plan['pins'].items():assert sha(p)==h,p
subprocess.run([sys.executable,'scripts/readback_domain_coverage_geometry.py','--plan',str(pp),'--output','metadata/domain_coverage_geometry_readback_20260927.json'],check=True)
