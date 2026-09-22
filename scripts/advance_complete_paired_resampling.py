#!/usr/bin/env python3
"""Gate full expanded paired branch resampling on audited point fits."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from assess_pae_sensitivity import checked_receipt
from readback_whole_proteome_family_coverage import sha


def require_values(actual, expected):
    for key,value in expected.items():
        if actual.get(key)!=value:
            raise ValueError('Completion binding mismatch: '+key)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Changed controller plan')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Changed dependency: '+path)
    verify()
    control=Path(plan['control']);output=Path(plan['resampling']);audit=Path(plan['audit'])
    if any(p.exists() for p in (control,output,audit)):raise FileExistsError('Use new immutable outputs')
    control.mkdir(parents=True)
    state=dict(status='waiting_for_exact_point_fit_controller',started_unix=time.time(),plan_sha256=plan_sha)
    def save():
        temporary=control/'state.tmp';temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(control/'state.json')
    save();dependency=plan['predecessor']
    while True:
        try:
            process=psutil.Process(dependency['pid'])
            live=process.create_time()==dependency['create_time'] and process.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify()
    inputs=Path(plan['inputs']);fits=Path(plan['fits']);fit_audit=Path(plan['fit_audit'])
    ir=checked_receipt(inputs);ar=checked_receipt(fit_audit);checked_receipt(Path(plan['models']))
    parent=json.loads(Path(plan['predecessor_receipt']).read_text())
    require_values(parent,dict(status='complete_expanded_supported_paired_point_fits_and_report_audit',
                   plan_sha256=sha(plan['predecessor_plan']),input_receipt_sha256=sha(inputs/'receipt.json'),
                   fit_receipt_sha256=sha(fits/'receipt.json'),audit_receipt_sha256=sha(fit_audit/'receipt.json'),
                   markers=plan['markers'],fits=4*plan['markers']))
    require_values(ar,dict(status='complete_paired_fit_audit',markers=plan['markers'],fits=4*plan['markers'],
                           fit_receipt_sha256=sha(fits/'receipt.json')))
    with (inputs/'marker_summary.tsv').open() as h:
        ready=[r for r in csv.DictReader(h,delimiter='\t') if r['status']=='ready_for_inference']
    if len(ready)!=plan['markers'] or ir['ready_markers']!=len(ready):raise ValueError('Incomplete marker scope')
    if min(int(r['retained_columns']) for r in ready)<max(plan['blocks']):raise ValueError('Block exceeds alignment length')
    draws=len(ready)*len(plan['blocks'])*plan['replicates']
    dimensions=dict(markers=len(ready),paired_draws=draws,maximum_fits=2*draws,
                    minimum_columns=min(int(r['retained_columns']) for r in ready),
                    resources=plan['resources'])
    (control/'actual_prelaunch_dimensions.json').write_text(json.dumps(dimensions,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    env['PATH']=str(Path(plan['iqtree']).parent)+os.pathsep+env.get('PATH','')
    commands=[('paired_resampling',[sys.executable,'scripts/run_paired_branch_resampling.py','--inputs',str(inputs),
               '--models',plan['models'],'--fits',str(fits),'--audit',str(fit_audit),'--output',str(output),
               '--replicates',str(plan['replicates']),'--blocks',*map(str,plan['blocks'])]),
              ('full_resampling_audit',[sys.executable,'scripts/summarize_paired_resampling.py','--inputs',str(inputs),
               '--fits',str(fits),'--resampling',str(output),'--output',str(audit)])]
    for stage,command in commands:
        verify()
        if shutil.disk_usage(control).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient free disk')
        if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient available memory')
        state.update(status='running_stage',stage=stage,command=command);save()
        with (control/(stage+'.log')).open('w') as log:
            subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    audited=checked_receipt(audit)
    require_values(audited,dict(status='complete_paired_resampling_audit',markers=len(ready),
                   batches=len(ready)*len(plan['blocks']),attempted_paired_draws=draws,
                   resampling_receipt_sha256=sha(output/'receipt.json')))
    if audited['validated_fits']!=2*(draws-audited['unestimable_draws']):raise ValueError('Incomplete resampling grid')
    verify()
    result=dict(status='complete_expanded_paired_branch_resampling_and_audit',plan_sha256=plan_sha,
                markers=len(ready),attempted_paired_draws=draws,validated_fits=audited['validated_fits'],
                unestimable_draws=audited['unestimable_draws'],fits_with_warnings=audited['fits_with_warnings'],
                predecessor_receipt_sha256=sha(plan['predecessor_receipt']),
                resampling_receipt_sha256=sha(output/'receipt.json'),audit_receipt_sha256=sha(audit/'receipt.json'),
                scope='Conditional paired site/block branch sampling sensitivity on fixed AA topologies. Not total calibrated uncertainty, independent biological coupling, selection or acceleration evidence.')
    (control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state.update(status=result['status']);save()


if __name__=='__main__':main()
