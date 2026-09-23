#!/usr/bin/env python3
"""Resume unchanged paired fits and audit after reboot; preserve source plan."""
import argparse
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


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    old_path=Path(plan['original_plan']);old=json.loads(old_path.read_text())
    if Path('/proc/sys/kernel/random/boot_id').read_text().strip()!=plan['recovery_boot_id']:raise ValueError('Recovery boot changed')
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Changed recovery plan')
        for path,digest in {**old['pins'],**plan['pins']}.items():
            if sha(path)!=digest:raise ValueError('Changed dependency: '+path)
    verify();control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
    original_control=Path(old['control']);fits=Path(old['fits']);audit=Path(old['audit'])
    if (original_control/'receipt.json').exists() or audit.exists():raise FileExistsError('Inspect completed or partial fit audit first')
    shutil.copy2(original_control/'state.json',control/'original_state.json')
    inputs=Path(old['inputs']);models=Path(old['models']);ir=checked_receipt(inputs);checked_receipt(models)
    readback=json.loads(Path(old['input_readback']).read_text())
    if readback['status']!='passed_complete_paired_inputs_from_qualified_arrays_readback' or readback['source_receipt_sha256']!=sha(inputs/'receipt.json'):raise ValueError('Input readback not bound')
    config=json.loads((fits/'config.json').read_text())
    if config['input_receipt_sha256']!=sha(inputs/'receipt.json') or config['model_receipt_sha256']!=sha(models/'receipt.json'):raise ValueError('Existing fit provenance differs')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    env['PATH']=str(Path(old['iqtree']).parent)+os.pathsep+env.get('PATH','')
    commands=[('checkpointed_fits',[sys.executable,'scripts/run_paired_marker_fits.py','--inputs',str(inputs.resolve()),'--models',str(models.resolve()),'--output',str(fits.resolve()),'--alrt','1000','--bootstrap','1000']),
              ('fit_audit',[sys.executable,'scripts/summarize_paired_marker_fits.py','--inputs',str(inputs.resolve()),'--models',str(models.resolve()),'--fits',str(fits.resolve()),'--output',str(audit.resolve())])]
    for stage,command in commands:
        verify()
        if psutil.virtual_memory().available<32*2**30 or shutil.disk_usage(control).free<100*2**30:raise ValueError('Insufficient resources')
        (control/'state.json').write_text(json.dumps(dict(status='running_stage',stage=stage,command=command,plan_sha256=plan_sha),indent=2)+'\n')
        with (control/(stage+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    fr=json.loads((fits/'receipt.json').read_text());ar=checked_receipt(audit);n=ir['ready_markers']
    if fr['status']!='complete_matched_topology_point_estimates' or fr['markers']!=n or ar['status']!='complete_paired_fit_audit' or ar['markers']!=n or ar['fits']!=4*n or ar['fit_receipt_sha256']!=sha(fits/'receipt.json'):raise ValueError('Incomplete recovered fit grid')
    verify()
    result=dict(status='complete_expanded_supported_paired_point_fits_and_report_audit',markers=n,fits=ar['fits'],paired_branches=ar['paired_branches'],fits_with_warnings=ar['fits_with_warnings'],plan_sha256=sha(old_path),input_receipt_sha256=sha(inputs/'receipt.json'),fit_receipt_sha256=sha(fits/'receipt.json'),audit_receipt_sha256=sha(audit/'receipt.json'),recovery_plan_sha256=plan_sha,scope='Original fit plan completed through checkpoint recovery; full point/report audit, not complete branch uncertainty or biological inference.')
    (original_control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    (control/'receipt.json').write_text(json.dumps(dict(status='complete_paired_fit_checkpoint_recovery',plan_sha256=plan_sha,original_controller_receipt_sha256=sha(original_control/'receipt.json')),indent=2)+'\n')


if __name__=='__main__':main()
