#!/usr/bin/env python3
"""Audit a complete BAli-Phy diagnostic batch only after verified producer exit."""
import argparse,json,subprocess,time
from pathlib import Path
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());digest=sha(a.plan)
    def verify():
        assert sha(a.plan)==digest
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();identity=json.loads(Path(plan['launch']).read_text())
    while True:
        verify()
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'-p','MainPID','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState']=='active' and int(state['MainPID'])==identity['pid'],state
        try:
            process=psutil.Process(identity['pid']);assert process.create_time()==identity['created'] and process.cmdline()==identity['command']
        except psutil.NoSuchProcess:
            time.sleep(1);continue
        time.sleep(30)
    verify();source_plan=json.loads(Path(plan['producer_plan']).read_text());root=Path(source_plan['output'])
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['jobs']==plan['expected_jobs'] and len(receipt['job_receipts'])==plan['expected_jobs'] and receipt['plan_sha256']==sha(plan['producer_plan'])
    for path,h in receipt['job_receipts'].items():assert sha(root/path)==h,path
    subprocess.run([plan['python'],plan['audit_script'],'--output',plan['output']],check=True)
    arp=Path(plan['output'])/'receipt.json';audit=json.loads(arp.read_text())
    assert audit['status']==plan['expected_audit_status'] and audit['expected_jobs']==plan['expected_jobs']
    assert sum(audit['dispositions'].values())==plan['expected_jobs'] and not audit['dispositions'].get('pending_receipt',0)
    for path,h in audit['artifacts'].items():assert sha(arp.parent/path)==h,path
    result=dict(status='complete_diagnostic_disposition_readback_review_outcomes',terminal_state=state,producer_receipt_sha256=sha(root/'receipt.json'),audit_receipt=str(arp),audit_receipt_sha256=sha(arp),dispositions=audit['dispositions'],handoff_plan_sha256=digest,scope='Complete accounting and declared diagnostic checks only. Failed jobs remain unresolved. No convergence, posterior qualification or full-project completion claim.')
    target=Path(plan['completion_metadata']);assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
