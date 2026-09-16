#!/usr/bin/env python3
"""Run the full discovery queue followed by independent native output readback."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import psutil
from run_clade_discovery import save,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());root=Path(plan['output']);root.mkdir(parents=True,exist_ok=True)
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Changed pinned input '+name)
    state=dict(status='running_discovery',pid=os.getpid(),created=psutil.Process().create_time(),plan_sha256=sha(a.plan))
    save(root/'controller_state.json',state)
    try:
        subprocess.run([plan['controller_python'],'-u','scripts/run_clade_discovery.py','--plan',str(a.plan)],check=True,cwd=plan['working_directory'])
        state['status']='running_independent_readback';save(root/'controller_state.json',state)
        subprocess.run([plan['orthofinder_python'],'scripts/readback_clade_discovery.py','--plan',str(a.plan),'--output',str(root/'readback.json')],check=True,cwd=plan['working_directory'],env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
        readback=json.loads((root/'readback.json').read_text())
        if readback['status']!='passed_complete_native_clade_discovery_readback' or readback['plan_sha256']!=sha(a.plan):raise ValueError('Readback did not pass')
        state.update(status='complete_discovery_and_readback',readback_sha256=sha(root/'readback.json'));save(root/'controller_state.json',state)
    except BaseException as exc:
        state.update(status='failed_requires_review',error=repr(exc));save(root/'controller_state.json',state);raise


if __name__=='__main__':main()
