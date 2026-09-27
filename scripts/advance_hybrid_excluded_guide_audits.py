#!/usr/bin/env python3
"""Wait for the full hybrid-excluded guide pair and audit every bootstrap tree."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());digest=sha(a.plan)
    root=Path(plan['controller_output']);root.mkdir(parents=True,exist_ok=False)
    def verify():
        assert sha(a.plan)==digest
        for path,value in plan['pins'].items():assert sha(path)==value,path
    verify();identity=plan['producer']
    while True:
        try:
            p=psutil.Process(identity['pid'])
            if p.create_time()!=identity['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==identity['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_two_full_hybrid_excluded_guides',identity['pid'],flush=True);time.sleep(30)
    verify();runs=Path(plan['runs']);r=json.loads((runs/'receipt.json').read_text())
    assert r['status']=='complete_two_hybrid_excluded_guide_executions_pending_audits'
    assert r['plan_sha256']==sha(plan['producer_plan'])
    assert len(r['runs'])==2 and {x['label'] for x in r['runs']}=={'profile','mafft'}
    proofs=[]
    for run in r['runs']:
        label=run['label'];folder=runs/label;output=Path(plan['audit_outputs'][label])
        assert sha(folder/'receipt.json')==run['receipt_sha256']
        command=[sys.executable,'scripts/audit_hybrid_excluded_species_guides.py','--run',str(folder),'--plan',plan['producer_plan'],'--output',str(output)]
        with (root/(label+'.log')).open('x') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
        verify();audit=json.loads((output/'receipt.json').read_text())
        assert audit['status']=='passed_full_hybrid_excluded_guide_bootstrap_readback'
        assert audit['taxa']==524 and audit['bootstrap_trees']==1000
        assert audit['source_receipt_sha256']==run['receipt_sha256']
        for name,value in audit['artifacts'].items():assert sha(output/name)==value
        proofs.append(dict(label=label,audit_receipt_sha256=sha(output/'receipt.json')))
    result=dict(status='complete_two_hybrid_excluded_guide_audits',runs=proofs,plan_sha256=digest,source_receipt_sha256=sha(runs/'receipt.json'),scope='Full saved tree and bootstrap-support readback; no fresh likelihood or SH-aLRT calculation, mixture-model adequacy or final species-tree claim.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
