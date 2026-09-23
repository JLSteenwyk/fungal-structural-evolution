#!/usr/bin/env python3
"""Run full tree-path readback after an exact completed benchmark controller."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from readback_tree_path_geometry_points import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Plan changed')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed dependency: '+path)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    if Path(plan['readback_output']).exists():raise FileExistsError('Existing readback output')
    state=out/'state.json';state.write_text(json.dumps({'status':'waiting_for_benchmark'})+'\n')
    dep=plan['predecessor']
    while True:
        try:
            p=psutil.Process(dep['pid']);live=p.create_time()==dep['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();cp=json.loads(Path(plan['benchmark_plan']).read_text());rp=Path(cp['output'])/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if (r['status']!='complete_all_cohort_descriptive_tree_path_geometry_benchmark'
            or r['plan_sha256']!=sha(plan['benchmark_plan']) or r['markers']!=plan['markers'] or r['pairs']!=plan['pairs']):raise ValueError('Incomplete or wrong benchmark')
    target=Path(cp['stages'][0]['output'])/'receipt.json'
    if r['stages'][0]['name']!='paths' or Path(r['stages'][0]['receipt'])!=target or sha(target)!=r['stages'][0]['sha256']:raise ValueError('Benchmark result binding differs')
    if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30 or shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient resources')
    state.write_text(json.dumps({'status':'checking_all_paths_and_inherited_fields'})+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    with (out/'readback.log').open('w') as log:subprocess.run([sys.executable,*plan['arguments']],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    rr=Path(plan['readback_output'])/'receipt.json';audit=json.loads(rr.read_text())
    if (audit['status']!='passed_all_tree_path_geometry_point_rows' or audit['markers']!=plan['markers']
            or sum(audit['counts'].values())!=plan['pairs'] or audit['tree_paths_checked']!=4*plan['pairs']
            or audit['source_hashes'][str(target)]!=sha(target) or sha(rp)!=rh):raise ValueError('Readback scope or source changed')
    verify();result=dict(status='complete_full_tree_path_geometry_readback_handoff',plan_sha256=ph,benchmark_controller_receipt_sha256=rh,readback_receipt_sha256=sha(rr),markers=audit['markers'],pairs=plan['pairs'],tree_paths_checked=audit['tree_paths_checked'])
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state.write_text(json.dumps({'status':result['status']})+'\n')


if __name__=='__main__':main()
