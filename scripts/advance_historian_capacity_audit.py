#!/usr/bin/env python3
"""Wait on the verified capacity producer, then audit all job dispositions."""
import json, subprocess, time
from pathlib import Path
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/historian_capacity_audit_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify()
    launch=json.loads(Path(plan['launch']).read_text())
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus','-p','MainPID'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid'],state
        try:
            process=psutil.Process(launch['pid'])
            assert process.create_time()==launch['create_time'] and process.cmdline()==launch['command']
        except psutil.NoSuchProcess:
            # Process may have exited between systemd readback and identity check.
            time.sleep(1);continue
        time.sleep(30)
    verify()
    capacity=json.loads(Path(plan['capacity_plan']).read_text());root=Path(capacity['output'])
    completed=json.loads((root/'receipt.json').read_text())
    assert completed['jobs']==324 and completed['plan_sha256']==sha(plan['capacity_plan'])
    assert len(completed['job_receipts'])==324
    for p,h in completed['job_receipts'].items():assert sha(root/p)==h,p
    subprocess.run([plan['python'],plan['audit_script'],'--output',plan['output']],check=True)
    audit=json.loads((Path(plan['output'])/'receipt.json').read_text())
    assert audit['status']=='all_job_dispositions_audited' and audit['expected_jobs']==324
    for p,h in audit['artifacts'].items():assert sha(Path(plan['output'])/p)==h
    receipt=dict(status='full_capacity_disposition_audit_complete_review_outcomes',producer_terminal_state=state,producer_receipt_sha256=sha(root/'receipt.json'),audit_receipt_sha256=sha(Path(plan['output'])/'receipt.json'),audit_dispositions=audit['dispositions'],plan_sha256=sha(pp),scope='All expected outcomes accounted for. Failed capacity outcomes remain unresolved; passing readback does not establish fitted ancestral reliability.')
    Path(plan['completion_metadata']).write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
