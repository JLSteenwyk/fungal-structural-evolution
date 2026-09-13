#!/usr/bin/env python3
"""Finish the frozen dataset acquisition, then launch full-ingroup lineage QC."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from advance_expanded_structural_confidence import process_identity
from prepare_lineage_busco_datasets import ROOT, sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download-pid',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new pipeline control output')
    identity=process_identity(a.download_pid)
    recovery=json.loads((ROOT/'metadata/lineage_busco_dataset_finalization_recovery.json').read_text())
    if recovery['original_pid']!=a.download_pid:raise ValueError('Downloader differs from recorded recovery')
    if identity is not None and 'scripts/prepare_lineage_busco_datasets.py' not in identity['command']:raise ValueError('Downloader PID was reused')
    a.output.mkdir(parents=True)
    scripts={name:sha(ROOT/'scripts'/name) for name in ['prepare_lineage_busco_datasets.py','run_lineage_busco_qc.py']}
    config={'download_pid':a.download_pid,'download_identity':identity,'script_sha256':sha(Path(__file__)),'stage_scripts_sha256':scripts,'resource_plan_sha256':sha(ROOT/'metadata/lineage_busco_qc_resource_plan.json'),'recovery_plan_sha256':sha(ROOT/'metadata/lineage_busco_dataset_finalization_recovery.json')}
    (a.output/'config.json').write_text(json.dumps(config,indent=2)+'\n');last=0
    while current:=process_identity(a.download_pid):
        if current!=identity:raise ValueError('Downloader PID identity changed')
        if time.monotonic()-last>=60:print('Waiting for original dataset download process',a.download_pid,flush=True);last=time.monotonic()
        time.sleep(10)
    stages=[('prepare_lineage_busco_datasets.py',['--output','data/busco-lineage-qc-v1','--resume'],'dataset_finalize.log'),('run_lineage_busco_qc.py',['--datasets','data/busco-lineage-qc-v1','--output','results/busco-lineage-v1'],'lineage_qc.log')]
    for script,args,log in stages:
        if sha(ROOT/'scripts'/script)!=scripts[script]:raise ValueError('Stage code changed')
        print('Starting',script,flush=True)
        with (a.output/log).open('w') as f:subprocess.run([sys.executable,str(ROOT/'scripts'/script),*args],cwd=ROOT,stdout=f,stderr=f,check=True)
    r=json.loads((ROOT/'results/busco-lineage-v1/receipt.json').read_text())
    if r['status']!='complete_full_ingroup_lineage_qc' or r['jobs']!=501 or r['completed']!=501:raise ValueError('Full lineage QC did not complete')
    result={'status':'complete_lineage_qc_pipeline','config_sha256':sha(a.output/'config.json'),'dataset_receipt_sha256':sha(ROOT/'data/busco-lineage-qc-v1/receipt.json'),'qc_receipt_sha256':sha(ROOT/'results/busco-lineage-v1/receipt.json'),'interpretation':'All 501 ingroup protein QC runs completed; comparative interpretation and annotation/contamination reviews remain separate.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
