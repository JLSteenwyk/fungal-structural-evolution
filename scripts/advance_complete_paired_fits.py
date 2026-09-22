#!/usr/bin/env python3
"""Run all eligible expanded paired fits after complete input readback."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from audit_busco_gene_copies import ROOT,sha,read_table
from assess_pae_sensitivity import checked_receipt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:raise ValueError('Plan changed')
        for name,digest in plan['pins'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Changed dependency: '+name)
    verify()
    control=ROOT/plan['control'];fits=ROOT/plan['fits'];audit=ROOT/plan['audit']
    if any(p.exists() for p in [control,fits,audit]):raise FileExistsError('Use new immutable output paths')
    control.mkdir(parents=True)
    state=dict(status='waiting_for_exact_paired_input_controller',started_unix=time.time(),plan_sha256=plan_sha)
    def save():
        p=control/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(control/'state.json')
    save()
    while True:
        try:
            proc=psutil.Process(plan['predecessor_pid'])
            live=proc.create_time()==plan['predecessor_create_time'] and proc.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify()
    previous=json.loads((ROOT/plan['predecessor_receipt']).read_text())
    inputs=ROOT/plan['inputs'];models=ROOT/plan['models']
    ir=checked_receipt(inputs);model_receipt=checked_receipt(models)
    readback_path=ROOT/plan['input_readback'];readback=json.loads(readback_path.read_text())
    if (previous['status']!='complete_full_cohort_paired_inputs_and_independent_array_readback'
            or previous['plan_sha256']!=sha(ROOT/plan['predecessor_plan'])
            or previous['input_receipt_sha256']!=sha(inputs/'receipt.json')
            or previous['readback_sha256']!=sha(readback_path)
            or readback['status']!='passed_complete_paired_inputs_from_qualified_arrays_readback'
            or readback['source_receipt_sha256']!=sha(inputs/'receipt.json')
            or model_receipt['status']!='complete_published_3di_model_validation'):
        raise ValueError('Verified paired inputs and published models required')
    ready=[r for r in read_table(inputs/'marker_summary.tsv') if r['status']=='ready_for_inference']
    if not ready or len(ready)!=ir['ready_markers'] or len(ready)>125:
        raise ValueError('Eligible marker scope differs')
    dimensions=dict(markers=len(ready),fits=4*len(ready),taxa_range=[min(int(r['eligible_taxa']) for r in ready),max(int(r['eligible_taxa']) for r in ready)],
                    columns_range=[min(int(r['retained_columns']) for r in ready),max(int(r['retained_columns']) for r in ready)],
                    input_receipt_sha256=sha(inputs/'receipt.json'),resources=plan['resources'])
    (control/'actual_prelaunch_dimensions.json').write_text(json.dumps(dimensions,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    env['PATH']=str(Path(plan['iqtree']).parent)+os.pathsep+env.get('PATH','')
    commands=[('supported_paired_fits',[sys.executable,'scripts/run_paired_marker_fits.py','--inputs',str(inputs),'--models',str(models),
               '--output',str(fits),'--alrt','1000','--bootstrap','1000']),
              ('fit_summary_audit',[sys.executable,'scripts/summarize_paired_marker_fits.py','--inputs',str(inputs),'--models',str(models),
                                    '--fits',str(fits),'--output',str(audit)])]
    for stage,command in commands:
        verify()
        if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
        if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient available memory')
        state.update(status='running_stage',stage=stage,command=command);save()
        with (control/(stage+'.log')).open('w') as f:
            subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    fr=json.loads((fits/'receipt.json').read_text());ar=checked_receipt(audit)
    if (fr['status']!='complete_matched_topology_point_estimates' or fr['markers']!=len(ready)
            or ar['status']!='complete_paired_fit_audit' or ar['markers']!=len(ready) or ar['fits']!=4*len(ready)
            or ar['fit_receipt_sha256']!=sha(fits/'receipt.json')):
        raise ValueError('Incomplete fit grid or output audit')
    verify()
    result=dict(status='complete_expanded_supported_paired_point_fits_and_report_audit',markers=len(ready),fits=ar['fits'],
                paired_branches=ar['paired_branches'],fits_with_warnings=ar['fits_with_warnings'],plan_sha256=plan_sha,
                input_receipt_sha256=sha(inputs/'receipt.json'),fit_receipt_sha256=sha(fits/'receipt.json'),audit_receipt_sha256=sha(audit/'receipt.json'),
                scope='All eligible markers fit with 1000 SH-aLRT and 1000 UFBoot sequence support, and three structural-model sensitivities on each AA topology. Output/report checks do not recompute likelihoods. Shared genealogy, branch uncertainty, direct geometry, prediction circularity and biological eligibility require further analyses.')
    (control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state.update(status=result['status']);save()


if __name__=='__main__':main()
