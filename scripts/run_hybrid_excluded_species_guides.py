#!/usr/bin/env python3
"""Infer both full hybrid-excluded guide trees serially from a pinned plan."""
import argparse
import fcntl
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    root=Path(plan['output']);root.mkdir(parents=True,exist_ok=True)
    lock=(root/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def verify():
        assert sha(args.plan)==ph
        for p,d in plan['pins'].items():assert sha(p)==d,p
    verify();completed=[]
    for label in ['profile','mafft']:
        verify()
        assert psutil.virtual_memory().available>=64*1024**3
        assert shutil.disk_usage(root).free>=100*1024**3
        folder=root/label;folder.mkdir(exist_ok=True)
        command=plan['commands'][label];config=folder/'config.json'
        expected=dict(command=command,plan_sha256=ph)
        if config.exists():assert json.loads(config.read_text())==expected
        else:config.write_text(json.dumps(expected,indent=2)+'\n')
        receipt_path=folder/'receipt.json'
        if receipt_path.exists():
            receipt=json.loads(receipt_path.read_text());assert receipt['config_sha256']==sha(config) and receipt['returncode']==0
            for n,d in receipt['artifacts'].items():assert sha(folder/n)==d
        else:
            start=time.monotonic()
            # An interrupted identical run uses IQ-TREE's existing checkpoint.
            with (folder/'controller.log').open('a') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
            verify()
            required=['guide.treefile','guide.contree','guide.ufboot','guide.iqtree','guide.log']
            assert all((folder/n).is_file() for n in required)
            receipt=dict(status='complete_hybrid_excluded_guide_execution_pending_full_audit',returncode=0,label=label,config_sha256=sha(config),elapsed_seconds_this_invocation=time.monotonic()-start,artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.name not in ['receipt.json','controller.log']},scope='Fresh full 524-taxon LG+F+G4 guide inference with 1000 SH-aLRT and UFB replicates. Full result audit and mixture-model sensitivity remain required.')
            receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
        completed.append(dict(label=label,receipt_sha256=sha(receipt_path)))
    verify();(root/'receipt.json').write_text(json.dumps(dict(status='complete_two_hybrid_excluded_guide_executions_pending_audits',plan_sha256=ph,runs=completed),indent=2)+'\n')


if __name__=='__main__':main()
