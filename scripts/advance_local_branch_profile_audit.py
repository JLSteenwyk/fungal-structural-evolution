#!/usr/bin/env python3
"""Validate the complete local branch-profile grid after its producer exits."""
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
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());digest=sha(a.plan)
    root=Path(plan['controller_output']);root.mkdir(parents=True,exist_ok=False)
    def verify():
        assert sha(a.plan)==digest
        for path,value in plan['pins'].items():assert sha(path)==value,path
    verify();identity=plan['producer']
    while True:
        try:
            p=psutil.Process(identity['pid'])
            if p.create_time()!=identity['create_time'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==identity['command'],'Producer command changed'
        except psutil.NoSuchProcess:break
        print('waiting_for_complete_local_profiles',p.pid,flush=True);time.sleep(30)
    verify()
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',identity['unit'],'-p','MainPID','-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state['MainPID']=='0' and state['ActiveState']=='inactive' and state['ExecMainStatus']=='0',state
    profiles=Path(plan['profiles']);r=json.loads((profiles/'receipt.json').read_text())
    assert r['status']=='complete_full_case_branch_parameter_profiles_pending_summary_audit'
    assert r['cases']==len(r['case_receipts'])==1632
    assert r['optimized_fits']==r['fresh_saved_fit_readbacks']==13056
    config=json.loads((profiles/'config.json').read_text())
    assert r['config_sha256']==sha(profiles/'config.json')
    assert config['script_sha256']==sha('scripts/profile_genus_longest_branch_parameter.py')
    assert config['plan_sha256']==sha(plan['producer_plan'])
    source_hash=sha(profiles/'receipt.json')
    command=[sys.executable,'scripts/audit_local_branch_parameter_profiles.py','--profiles',str(profiles),'--fits',plan['fits'],'--slices',plan['slices'],'--plan',plan['producer_plan'],'--exposure',plan['exposure'],'--output',plan['audit_output']]
    with (root/'audit.log').open('x') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
    verify();out=Path(plan['audit_output']);audit=json.loads((out/'receipt.json').read_text())
    assert audit['status']=='passed_full_branch_parameter_profile_artifact_and_numeric_audit'
    assert audit['cases']==1632 and audit['optimized_points']==13056
    assert audit['source_receipt_sha256']==source_hash==sha(profiles/'receipt.json')
    for name,value in audit['artifacts'].items():assert sha(out/name)==value
    result=dict(status='complete_full_local_branch_profile_audit_handoff',cases=1632,points=13056,source_receipt_sha256=source_hash,audit_receipt_sha256=sha(out/'receipt.json'),plan_sha256=digest,producer_terminal_state=state,scope='All finite-grid saved parameters, likelihoods and artifacts checked. No global optimum, confidence interval or selection conclusion established.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
