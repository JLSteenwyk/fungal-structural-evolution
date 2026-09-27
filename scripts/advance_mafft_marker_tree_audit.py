#!/usr/bin/env python3
"""Audit the complete MAFFT marker batch after its identified producer exits."""
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
        print('waiting_for_full_mafft_marker_batch',identity['pid'],flush=True);time.sleep(30)
    verify();trees=Path(plan['trees']);source=trees/'receipt.json';receipt=json.loads(source.read_text())
    assert receipt['status']=='complete_mafft_marker_tree_execution_pending_full_audit'
    assert receipt['plan_sha256']==sha(plan['producer_plan'])
    assert receipt['markers']==len(receipt['results'])==125
    assert len({r['marker'] for r in receipt['results']})==125
    assert all(r['status']=='inferred' and r['returncode']==0 for r in receipt['results'])
    source_hash=sha(source)
    commands=[
      [sys.executable,'scripts/audit_mafft_marker_tree_support.py','--trees',str(trees),'--matrix',plan['matrix'],'--output',plan['support_output']],
      [sys.executable,'scripts/readback_marker_tree_snapshot.py','--trees',str(trees),'--matrix',plan['matrix'],'--snapshot',plan['support_output'],'--output',plan['readback_output']]]
    for index,command in enumerate(commands):
        verify()
        with (root/f'stage_{index}.log').open('x') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
    verify();support=Path(plan['support_output']);readback=Path(plan['readback_output'])
    sr=json.loads((support/'receipt.json').read_text());ar=json.loads((readback/'receipt.json').read_text())
    assert sr['status']=='complete_support_audit' and sr['completed_markers']==sr['planned_markers']==125 and not sr['pending_markers']
    assert ar['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and ar['markers']==125
    assert ar['snapshot_receipt_sha256']==sha(support/'receipt.json')
    for folder,r in [(support,sr),(readback,ar)]:
        for n,d in r['artifacts'].items():assert sha(folder/n)==d
    assert sha(source)==source_hash
    result=dict(status='complete_full_mafft_marker_tree_audit_handoff',markers=125,source_receipt_sha256=source_hash,support_receipt_sha256=sha(support/'receipt.json'),readback_receipt_sha256=sha(readback/'receipt.json'),plan_sha256=digest,internal_splits_checked=ar['internal_splits_checked'],retained_alignment_characters_checked=ar['retained_alignment_characters_checked'],scope='Every completed marker input reconstructed from source matrix and coverage rule; full supported split grid independently checked by graph-edge removal. No biological adequacy, preferred topology or cause of discordance established.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
