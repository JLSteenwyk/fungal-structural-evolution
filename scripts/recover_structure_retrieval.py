#!/usr/bin/env python3
"""Resume the existing AFDB retrieval runner with a free-disk guard."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from readback_whole_proteome_family_coverage import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Changed dependency: '+path)
    audit=json.loads(Path(plan['log_readback']).read_text())
    if sha(audit['path'])!=audit['sha256']:raise ValueError('Retrieval log changed before recovery')
    output=Path(plan['output']);output.mkdir(parents=True,exist_ok=False)
    if shutil.disk_usage(output).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk headroom')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    command=[sys.executable,'scripts/retrieve_matched_models.py'];stopped=False
    with (output/'retrieval.log').open('w') as handle:
        process=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        while process.poll() is None:
            if shutil.disk_usage(output).free<plan['resources']['emergency_free_disk_gib']*2**30:
                stopped=True;os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=180)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                break
            time.sleep(20)
    if process.returncode:raise RuntimeError('Retrieval exited '+str(process.returncode)+'; inspect preserved log')
    with (output/'retrieval.log').open() as handle:
        last=''
        for line in handle:
            if line.strip():last=line.strip()
    complete=last.startswith('Snapshot queue completed') and not stopped
    result={'status':'complete_snapshot_retrieval_queue' if complete else 'retrieval_stopped_pending_review',
            'plan_sha256':plan_sha,'stopped_for_disk_headroom':stopped,'last_log_line':last,
            'scope':'Retrieval runner completion only; individual records may contain errors or no-match outcomes. Not complete atlas coverage or coordinate-confidence qualification.'}
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
