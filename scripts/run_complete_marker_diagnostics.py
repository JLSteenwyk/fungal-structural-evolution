#!/usr/bin/env python3
"""Run complete marker-tree audits and support-sensitive guide diagnostics."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',required=True,type=Path)
    a=p.parse_args();plan=json.loads(a.plan.read_text());plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:raise ValueError('Plan changed')
        for name,digest in plan['pins'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Changed source: '+name)
    verify()
    out=ROOT/plan['control']
    if out.exists() or any((ROOT/s['output']).exists() for s in plan['stages']):raise FileExistsError('Use fresh immutable diagnostics')
    out.mkdir(parents=True)
    state=dict(status='starting',plan_sha256=plan_sha,started_unix=time.time(),completed_stages=[])
    def save():
        p=out/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(out/'state.json')
    save()
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    for stage in plan['stages']:
        verify()
        if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
        command=[sys.executable]+stage['arguments']
        state.update(status='running_stage',stage=stage['name'],command=command);save()
        with (out/(stage['name']+'.log')).open('w') as f:
            subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        receipt_path=ROOT/stage['receipt'];result=json.loads(receipt_path.read_text())
        if any(result.get(k)!=v for k,v in stage['expected'].items()):raise ValueError('Unexpected stage status/scope')
        state['completed_stages'].append(dict(name=stage['name'],receipt=stage['receipt'],sha256=sha(receipt_path)));save()
    verify()
    result=dict(status='complete_full_marker_support_and_guide_diagnostics',markers=125,plan_sha256=plan_sha,stages=state['completed_stages'],
                scope='All planned supported marker trees audited; full-marker homogeneous-guide conflict and role-separation diagnostics. SH-aLRT is not bootstrap support; restricted guide rows are dependent. Does not establish final species topology, rooting, reconciliation or biological causes of discordance.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state.update(status=result['status']);save()


if __name__=='__main__':main()
