#!/usr/bin/env python3
"""Execute pinned collection site-rate exports and full output audits."""
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
verify()
verified={}
for expected in plan['audits']:
    rp=Path(expected['path'])/'receipt.json';r=json.loads(rp.read_text())
    if any(r[k]!=v for k,v in expected['counts'].items()):raise ValueError('Unexpected rate audit scope')
    for name,h in r['artifacts'].items():
        if sha(rp.parent/name)!=h:raise ValueError('Changed audited rate artifact')
    verified[str(rp)]=sha(rp)
(control/'receipt.json').write_text(json.dumps(dict(status='complete_collection_site_rate_exports_and_audits',plan_sha256=ph,audit_receipts=verified),indent=2)+'\n')
