#!/usr/bin/env python3
"""Wait for exact numeric-readback processes, then summarize each complete workload."""
import argparse,json,subprocess,sys,time
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
def verify():
    if sha(a.plan)!=ph:raise ValueError('Changed plan')
    for path,h in plan['pins'].items():
        if sha(path)!=h:raise ValueError('Changed pin '+path)
verify();control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
for stage in plan['stages']:
    dep=stage['dependency'];(control/'state.json').write_text(json.dumps(dict(stage=stage['name'],status='waiting'))+'\n')
    while psutil.pid_exists(dep['pid']):
        try:
            proc=psutil.Process(dep['pid'])
            if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
            if proc.cmdline()!=dep['cmdline']:raise ValueError('Dependency identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify()
    audit=json.loads((Path(stage['readback'])/'receipt.json').read_text())
    if audit['plan_sha256']!=sha(stage['readback_plan']):raise ValueError('Readback plan differs')
    command=[sys.executable,'scripts/summarize_duplication_alignment_orders.py']
    for key in ['producer','readback','output']:command.extend(['--'+key,stage[key]])
    with (control/(stage['name']+'.log')).open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
verify();(control/'receipt.json').write_text(json.dumps(dict(status='complete_primary_reference_order_summaries',plan_sha256=ph,receipts={s['output']:sha(Path(s['output'])/'receipt.json') for s in plan['stages']}),indent=2)+'\n')
