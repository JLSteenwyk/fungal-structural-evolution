#!/usr/bin/env python3
"""Join completed full-cohort geometry to audited matched-topology tree paths."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from readback_domain_boundary_clusters import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',required=True,type=Path)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan);pins=dict(plan['pins']);pins[str(a.plan)]=ph
    def verify():
        for p,h in pins.items():
            if sha(p)!=h:raise ValueError('Changed dependency: '+p)
    verify();out=Path(plan['output'])
    if out.exists() or any(Path(s['output']).exists() for s in plan['stages']):raise FileExistsError('Use fresh output paths')
    out.mkdir(parents=True);completed=[]
    def state(status,**extra):
        p=out/'state.tmp';p.write_text(json.dumps(dict(status=status,completed_stages=completed,**extra),indent=2)+'\n');p.replace(out/'state.json')
    def receipt(folder):
        path=Path(folder)/'receipt.json';pins[str(path)]=sha(path);return json.loads(path.read_text())
    try:
        state('waiting_for_geometry_and_paired_fit_controllers')
        for dependency in plan['dependencies']:
            while True:
                try:
                    p=psutil.Process(dependency['pid']);live=p.create_time()==dependency['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
                except psutil.NoSuchProcess:live=False
                if not live:break
                time.sleep(20)
            r=receipt(dependency['output'])
            if r['status']!=dependency['status'] or r['plan_sha256']!=sha(dependency['plan']):raise ValueError('Unsuccessful or wrong predecessor')
        verify();gp=json.loads(Path(plan['geometry_plan']).read_text());fp=json.loads(Path(plan['fit_plan']).read_text())
        gc=receipt(gp['output']);fc=receipt(fp['control']);fit=receipt(fp['fits']);audit=receipt(fp['audit'])
        geometry=receipt(plan['geometry']);ga=receipt(plan['geometry_audit']);inputs_sha=sha(Path(plan['inputs'])/'receipt.json')
        if (fc['status']!='complete_expanded_supported_paired_point_fits_and_report_audit'
                or fc['plan_sha256']!=sha(plan['fit_plan']) or fc['input_receipt_sha256']!=inputs_sha
                or fc['fit_receipt_sha256']!=sha(Path(fp['fits'])/'receipt.json')
                or fc['audit_receipt_sha256']!=sha(Path(fp['audit'])/'receipt.json')
                or fit['status']!='complete_matched_topology_point_estimates'
                or fit['markers']!=plan['markers'] or audit['markers']!=plan['markers']
                or audit['status']!='complete_paired_fit_audit' or audit['fits']!=4*plan['markers']
                or audit['fit_receipt_sha256']!=sha(Path(fp['fits'])/'receipt.json')
                or geometry['status']!='complete_paired_site_geometry' or geometry['markers']!=plan['markers']
                or geometry['source_receipts']['inputs']['sha256']!=inputs_sha
                or ga['status']!='passed_complete_paired_grid_character_and_sampled_geometry_readback'
                or ga['pair_rows']!=plan['pairs'] or ga['source_receipts']['inputs']!=inputs_sha
                or ga['source_receipts']['comparisons']!=sha(Path(plan['geometry'])/'receipt.json')):
            raise ValueError('Incomplete fit/geometry scope or wrong source binding')
        for stage in gc['stages']:
            if sha(stage['receipt'])!=stage['sha256']:raise ValueError('Geometry completion artifact changed')
        recovery=json.loads(Path(plan['dependencies'][1]['output'],'receipt.json').read_text())
        if recovery['original_controller_receipt_sha256']!=sha(Path(fp['control'])/'receipt.json'):raise ValueError('Wrong recovered fit controller')
        if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient memory')
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
        for stage in plan['stages']:
            verify()
            if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
            command=[sys.executable,*stage['arguments']];state('running_stage',stage=stage['name'])
            with (out/(stage['name']+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            r=receipt(stage['output'])
            if r['status']!=stage['status']:raise ValueError('Incomplete benchmark stage')
            if stage['name']=='paths' and (r['markers']!=plan['markers'] or r['accepted_pairs']!=ga['counts']['accepted'] or r['excluded_pairs']!=ga['counts']['excluded']):raise ValueError('Incomplete benchmark grid')
            completed.append(dict(name=stage['name'],receipt=str(Path(stage['output'])/'receipt.json'),sha256=sha(Path(stage['output'])/'receipt.json')))
        verify();result=dict(status='complete_all_cohort_descriptive_tree_path_geometry_benchmark',markers=plan['markers'],pairs=plan['pairs'],plan_sha256=ph,stages=completed,scope='Full matched-topology path/geometry point-estimate join and within-marker descriptive ranks. Numerical path checks sampled by original producer. Not independent pairwise tests, calibrated physical branch lengths, uncertainty analysis or biological acceleration. Figure requires visual review; full output readback remains required.')
        (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state(result['status'])
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
