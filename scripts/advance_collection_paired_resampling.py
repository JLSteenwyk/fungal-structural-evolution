#!/usr/bin/env python3
"""Execute a pinned full-cohort resampling plan followed by its artifact audit."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
import psutil
from paired_collection_resampling_sources import sha

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
plan=json.loads(a.plan.read_text());ph=sha(a.plan)
def verify():
    if sha(a.plan)!=ph:raise ValueError('Changed plan')
    for path,h in plan['pins'].items():
        if sha(Path(path))!=h:raise ValueError('Changed pinned source: '+path)
verify()
if shutil.disk_usage('.').free < plan['resources']['minimum_free_disk_gib']*2**30 or psutil.virtual_memory().available < plan['resources']['minimum_available_memory_gib']*2**30:raise RuntimeError('Insufficient resources')
control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
for stage in plan['stages']:
    verify();(control/'state.json').write_text(json.dumps(stage,indent=2)+'\n')
    with (control/(stage['name']+'.log')).open('w') as log:
        subprocess.run([sys.executable,*stage['arguments']],stdout=log,stderr=subprocess.STDOUT,check=True)
verify();r=json.loads((Path(plan['audit'])/'receipt.json').read_text())
if r['status']!='complete_paired_resampling_audit' or r['markers']!=125 or r['attempted_paired_draws']!=75000:raise ValueError('Unexpected audit scope')
(control/'receipt.json').write_text(json.dumps(dict(status='complete_collection_paired_resampling_and_audit',plan_sha256=ph,audit_receipt_sha256=sha(Path(plan['audit'])/'receipt.json')),indent=2)+'\n')
