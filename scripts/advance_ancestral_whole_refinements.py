#!/usr/bin/env python3
"""Wait for the full independent multistart audit before selecting refinement starts."""
import json,subprocess,time
from pathlib import Path
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/ancestral_whole_refinement_handoff_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();identity=json.loads(Path(plan['audit_launch']).read_text())
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',identity['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus','-p','MainPID'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState']=='active' and int(state['MainPID'])==identity['pid'],state
        try:
            p=psutil.Process(identity['pid']);assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline']
        except psutil.NoSuchProcess:
            time.sleep(1);continue
        time.sleep(30)
    verify();rp=Path(plan['audit_receipt']);r=json.loads(rp.read_text())
    assert r['status']=='complete_468_multistart_reports_and_independent_likelihoods' and r['fits']==468
    assert sha(r['source_receipt_path'])==r['source_receipt_sha256']
    for f,h in r['artifacts'].items():assert sha(rp.parent/f)==h,f
    closure=dict(r);closure.update(completed_receipt_path=str(rp),completed_receipt_sha256=sha(rp),verified_audit_terminal_state=state)
    cp=Path('metadata/ancestral_whole_multistart_audit_completed_20260927.json');assert not cp.exists();cp.write_text(json.dumps(closure,indent=2)+'\n')
    subprocess.run([plan['python'],'scripts/prepare_ancestral_whole_refinements.py'],check=True)
    prepared=json.loads(Path('metadata/ancestral_whole_refinement_plan_20260927.json').read_text());assert len(prepared['jobs'])==156
    subprocess.run([plan['python'],'scripts/run_ancestral_whole_refinements.py','--plan','metadata/ancestral_whole_refinement_plan_20260927.json'],check=True)
    print('whole_refinements_finished_pending_independent_audit',flush=True)
if __name__=='__main__':main()
