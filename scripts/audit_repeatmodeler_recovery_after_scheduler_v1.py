#!/usr/bin/env python3
"""Audit RepeatModeler recovery output only after its submitted Slurm array terminates."""
import argparse, hashlib, json, re, subprocess, time
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def terminal(states):
    allowed={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
    return bool(states) and all(x.split()[0].split('+')[0] in allowed for x in states)

def run(command): return subprocess.run(command,text=True,capture_output=True,check=False)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path);p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:
        assert terminal(['COMPLETED','FAILED+']) and not terminal([]) and not terminal(['PENDING'])
        print(json.dumps(dict(status='passed_repeatmodeler_recovery_audit_waiter_literal_control'),indent=2));return
    if a.plan is None:p.error('--plan is required unless --self-test is used')
    plan=json.loads(a.plan.read_text())
    for path,expected in plan['pins'].items():
        if sha(path)!=expected:raise ValueError('Changed pinned input: '+path)
    output=Path(plan['audit_receipt'])
    if output.exists():raise FileExistsError('Fresh audit receipt required')
    submission=Path(plan['submission_receipt'])
    while not submission.exists():
        print('waiting_for_repeatmodeler_recovery_submission_receipt',flush=True);time.sleep(plan['poll_seconds'])
    record=json.loads(submission.read_text())
    if record.get('status')!='submitted_repeatmodeler_recovery_array_after_verified_scheduler_return':raise ValueError('Unexpected submission receipt')
    for name in ('task_file','sbatch'):
        if record.get(name)!=plan[name] or record.get(name+'_sha256')!=plan[name+'_sha256']:raise ValueError('Submission input differs: '+name)
    match=re.fullmatch(r'Submitted batch job (\d+)',record.get('submission_stdout',''))
    if not match:raise ValueError('Cannot identify submitted recovery array')
    job=match.group(1)
    while True:
        ping=run(['scontrol','ping'])
        active=run(['squeue','--noheader','--jobs',job,'--format=%i']) if not ping.returncode else None
        if ping.returncode or active.returncode or active.stdout.strip():
            print('waiting_for_repeatmodeler_recovery_array_terminal',flush=True);time.sleep(plan['poll_seconds']);continue
        history=run(['sacct','--noheader','--parsable2','--jobs',job,'--format=JobIDRaw,State'])
        states=[line.split('|',1)[1] for line in history.stdout.splitlines() if '|' in line and line.split('|',1)[0].startswith(job)]
        if history.returncode or not terminal(states):
            print('waiting_for_repeatmodeler_recovery_array_accounting',flush=True);time.sleep(plan['poll_seconds']);continue
        result=run(['/home/bizon/anaconda3/bin/python',plan['audit_script'],'--config',plan['config'],'--receipt',str(output)])
        if result.returncode:raise RuntimeError('Full-panel audit failed: '+(result.stdout+result.stderr).strip())
        print(json.dumps(dict(status='completed_repeatmodeler_recovery_terminal_audit',job_id=job,terminal_states=states,audit_receipt=str(output)),indent=2));return
if __name__=='__main__':main()
