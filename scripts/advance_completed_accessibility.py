#!/usr/bin/env python3
"""Advance full-cohort ASA union, projection and normalization after complete source audits."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import psutil
from catalog_whole_proteome_structures import sha


def require_fields(result,expected):
    if any(result.get(k)!=v for k,v in expected.items()):raise ValueError('Incomplete or mismatched stage receipt')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed dependency: '+p)
    verify();out=Path(plan['control']);out.mkdir(parents=True,exist_ok=False)
    def state(s,**kw):(out/'state.json').write_text(json.dumps(dict(status=s,**kw),indent=2)+'\n')
    state('waiting_for_all_three_source_audits')
    while True:
        live=False
        for previous in plan['predecessors']:
            try:
                p=psutil.Process(previous['pid']);live|=p.create_time()==previous['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:pass
        if not live:break
        time.sleep(20)
    verify()
    for previous in plan['predecessors']:
        r=json.loads(Path(previous['receipt']).read_text())
        require_fields(r,{'status':'complete_full_accessibility_audit','config_sha256':sha(previous['config'])})
    if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30 or psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient resources')
    completed=[]
    for stage in plan['stages']:
        verify();state('running',stage=stage['name'],completed=completed)
        with (out/(stage['name']+'.log')).open('w') as log:
            subprocess.run(stage['command'],check=True,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
        rp=Path(stage['receipt']);r=json.loads(rp.read_text());require_fields(r,stage['expected'])
        completed.append({'name':stage['name'],'receipt':str(rp),'sha256':sha(rp)})
    verify()
    for r in completed:
        if sha(r['receipt'])!=r['sha256']:raise ValueError('Completed artifact changed')
    result={'status':'complete_all_cohort_accessibility_projection_and_normalization_readbacks','plan_sha256':ph,'completed_stages':completed,'models':plan['models'],'observed_sites':plan['observed_sites'],'scope':'All source audits, disjoint full-model ASA union, paired projection and normalization full-row readbacks completed. ASA method and prediction accuracy not independently revalidated here; controlled evolutionary and functional analyses remain separate.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state('complete')


if __name__=='__main__':main()
