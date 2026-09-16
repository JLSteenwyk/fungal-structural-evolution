#!/usr/bin/env python3
"""Wait for exact discovery completion, then merge and independently read back."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import psutil


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def save(path,data):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);a=ap.parse_args()
    config=json.loads(a.config.read_text());out=Path(config['controller_output']);out.mkdir(parents=True,exist_ok=False)
    state=dict(status='waiting_for_discovery_and_readback',pid=os.getpid(),created=psutil.Process().create_time(),config_sha256=sha(a.config));save(out/'state.json',state)
    def verify():
        for name,h in config['pins'].items():
            if sha(name)!=h:raise ValueError('Changed input '+name)
    try:
        verify();producer=config['producer']
        while True:
            try:
                p=psutil.Process(producer['pid'])
                active=abs(p.create_time()-producer['created'])<.05 and p.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:active=False
            if not active:break
            time.sleep(15)
        verify();discovery=Path(config['discovery_output']);ds=json.loads((discovery/'controller_state.json').read_text())
        if ds['status']!='complete_discovery_and_readback' or ds['plan_sha256']!=config['discovery_plan_sha256'] or ds['readback_sha256']!=sha(discovery/'readback.json'):
            raise ValueError('Successful complete discovery handoff required')
        if psutil.virtual_memory().available<config['minimum_available_memory_gib']*2**30 or shutil.disk_usage(out).free<config['minimum_free_disk_gib']*2**30:raise RuntimeError('Merge resource gate')
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        commands=[('merging',[config['controller_python'],'scripts/merge_clade_discovery_partitions.py','--plan',config['merge_plan']]),('independent_readback',[config['native_python'],'scripts/readback_guide_discovery_merge.py','--plan',config['merge_plan'],'--output',config['readback_output']])]
        for label,command in commands:
            state.update(status=label,command=command);save(out/'state.json',state)
            with (out/(label+'.log')).open('w') as log:
                subprocess.run(command,check=True,env=env,stdout=log,stderr=subprocess.STDOUT,cwd=config['working_directory'])
        result=json.loads(Path(config['readback_output']).read_text())
        if result['status']!='passed_complete_guide_discovery_merge_readback' or result['plan_sha256']!=sha(config['merge_plan']):raise ValueError('Independent merge readback did not pass')
        state.update(status='complete_merged_partitions_and_readback',readback_sha256=sha(config['readback_output']));save(out/'state.json',state);save(out/'receipt.json',state)
    except BaseException as exc:
        state.update(status='failed_requires_review',error=repr(exc));save(out/'state.json',state);raise


if __name__=='__main__':main()
