#!/usr/bin/env python3
"""Audit the fourth crossed PMSF run after its identified controller finishes."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());digest=sha(args.plan)
    root=Path(plan['controller_output']);root.mkdir(parents=True,exist_ok=False)
    lock=(root/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def verify():
        assert sha(args.plan)==digest
        for p,d in plan['pins'].items():assert sha(p)==d,p
    verify();identity=plan['producer']
    while True:
        try:
            p=psutil.Process(identity['pid'])
            if p.create_time()!=identity['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==identity['cmdline'],'Producer identity changed'
        except psutil.NoSuchProcess:break
        print('waiting_for_identified_fourth_pmsf_controller',identity['pid'],flush=True);time.sleep(30)
    verify();run=Path(plan['run']);receipt=json.loads((run/'receipt.json').read_text())
    assert receipt['status']=='complete_pmsf_execution_pending_full_audit' and receipt['returncode']==0
    source_hash=sha(run/'receipt.json')
    assert receipt['config_sha256']==sha(run/'config.json')
    command=[sys.executable,'scripts/audit_species_pmsf.py','--run',str(run),'--output',plan['audit_output']]
    with (root/'audit.log').open('x') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
    verify();out=Path(plan['audit_output']);proof=json.loads((out/'receipt.json').read_text())
    assert proof['status']=='passed_pmsf_profile_tree_and_bootstrap_readback'
    assert proof['source_receipt_sha256']==source_hash==sha(run/'receipt.json')
    assert proof['taxa']==526 and proof['sites']==63750 and proof['bootstrap_trees']==1000
    for name,d in proof['artifacts'].items():assert sha(out/name)==d
    result=dict(status='complete_fourth_pmsf_full_readback_handoff',plan_sha256=digest,source_receipt_sha256=source_hash,audit_receipt_sha256=sha(out/'receipt.json'),taxa=526,sites=63750,bootstrap_trees=1000,command=command,scope='Full saved profile/tree/bootstrap readback; all-four topology comparisons, rooting, adequacy, hybrid and other taxon/marker sensitivities remain required. No final species framework claim.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
