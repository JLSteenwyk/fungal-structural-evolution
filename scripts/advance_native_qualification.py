#!/usr/bin/env python3
"""Qualify a completed coordinate audit after its exact local PAE exporter finishes."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from advance_prediction_snapshot import ROOT, identity, sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True)
    a=p.parse_args();c=json.loads(a.config.read_text());config_sha=sha(a.config)
    control=ROOT/c['control_output'];control.mkdir(parents=True,exist_ok=True)
    lock=(control/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (control/'launch.json').exists() or (control/'receipt.json').exists():
        raise FileExistsError('Inspect existing execution before recovery')
    def verify():
        if sha(a.config)!=config_sha:raise ValueError('Changed configuration')
        for name,h in c['pinned_files'].items():
            if sha(ROOT/name)!=h:raise ValueError('Changed pinned dependency: '+name)
    verify();print('Waiting for local PAE exporter',c['predecessor_pid'],flush=True)
    while identity(c['predecessor_pid'])==c['predecessor_start_ticks']:time.sleep(20)
    verify();pae=ROOT/c['pae'];coordinates=ROOT/c['coordinates']
    pr=json.loads((pae/'receipt.json').read_text());cr=json.loads((coordinates/'receipt.json').read_text())
    if (pr['status']!='complete_local_mapping_bound_pae_export' or pr['models_verified']!=c['models']
            or pr['models_requested']!=c['models'] or pr['models_failed']!=0
            or cr['status']!='complete_native_3di_coordinate_audit' or cr['models']!=c['models']
            or pr['mapping_receipt_sha256']!=cr['mapping_receipt_sha256']):
        raise ValueError('Completed PAE and coordinate model/mapping universes differ')
    if (ROOT/c['output']).exists():raise FileExistsError('Use a fresh qualification output')
    if shutil.disk_usage(ROOT).free<c['resource_plan']['output_allowance_bytes']:raise RuntimeError('Insufficient output headroom')
    command=[sys.executable,'scripts/qualify_native_pae.py','--coordinates',c['coordinates'],
             '--pae',c['pae'],'--output',c['output']]
    launch={'command':command,'config_sha256':config_sha,'pae_receipt_sha256':sha(pae/'receipt.json'),
            'coordinate_receipt_sha256':sha(coordinates/'receipt.json'),'resource_plan':c['resource_plan']}
    (control/'launch.json').write_text(json.dumps(launch,indent=2)+'\n')
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    with (control/'qualification.log').open('w') as log:
        subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    verify();result={'status':'completed_native_qualification_command_pending_readback',
        'config_sha256':config_sha,'qualification_receipt_sha256':sha(ROOT/c['output']/'receipt.json'),
        'interpretation':'Qualification command completed after source checks; independent output/context readback and paired-input preparation remain separate.'}
    (control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
