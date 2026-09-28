#!/usr/bin/env python3
"""Queue every four-chain categorical report after producer/extractor handoff."""
import argparse
from collections import defaultdict,Counter
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import time
import psutil
from ancestral_chain_attempt import run_attempt,sha,write_json


def grouped_jobs(jobs):
    groups=defaultdict(list)
    for job in jobs:groups[job['config']['model_input_identity']].append(job)
    for group,rows in groups.items():
        assert len(rows)==len({j['chain']['chain_id'] for j in rows})==len({j['chain']['seed'] for j in rows})==4,group
    return dict(groups)


def ready(rows,source):
    for job in rows:
        path=Path(source)/job['chain']['chain_id']/'disposition.json'
        if not path.exists():return False
        if json.loads(path.read_text())['status']!='state_trace_complete_not_posterior_qualification':return False
    return True


def dependency_state(launch):
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    if state['ActiveState'] in ['inactive','failed']:return state
    assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
    try:
        p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
    except psutil.NoSuchProcess:
        return None  # Re-poll transitional disappearance; never restart producer.
    return state


def process_group(group,plan):
    target=Path(plan['output']).resolve()/group;target.mkdir(parents=True,exist_ok=True)
    pins={str(Path(p).absolute()):h for p,h in plan['pins'].items()}
    command=['/usr/bin/env','-C',plan['working_directory'],'/usr/bin/prlimit','--as='+str(plan['per_process_address_space_bytes']),'--',
        plan['python'],str(Path('scripts/prepare_ancestral_state_quartet_report.py').resolve()),
        '--producer-plan',str(Path(plan['producer_plan']).resolve()),'--extraction-plan',str(Path(plan['extraction_plan']).resolve()),
        '--group',group,'--output','{attempt}/report','--python',str(Path(plan['diagnostic_python']).absolute())]
    config=dict(command=command,timeout_seconds=plan['per_group_timeout_seconds'],pins=pins,group=group)
    attempt=run_attempt(target/'processing',config);ar=json.loads(attempt.read_text())
    result=dict(group=group,attempt_receipt=str(attempt),attempt_receipt_sha256=sha(attempt))
    if ar['exit_code']!=0:
        result['status']='categorical_report_attempt_failed_requires_review'
    else:
        receipt=attempt.parent/'report/receipt.json';r=json.loads(receipt.read_text())
        assert r['status']=='verified_quartet_categorical_reports_complete_not_posterior_qualification' and r['group']==group
        assert sha(r['report'])==r['report_sha256'] and sha(r['manifest'])==r['manifest_sha256']
        rr=json.loads(Path(r['report']).read_text())
        assert set(rr['outputs'])=={'250','500'}
        for name,h in rr['artifacts'].items():assert sha(Path(r['report']).parent/name)==h
        result.update(status='categorical_reports_complete_not_posterior_qualification',receipt=str(receipt),receipt_sha256=sha(receipt))
    write_json(target/'disposition.json',result)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify()
    producer=json.loads(Path(plan['producer_plan']).read_text());extraction=json.loads(Path(plan['extraction_plan']).read_text())
    launches=[json.loads(Path(p).read_text()) for p in plan['dependency_launches']]
    for launch,path in zip(launches,[plan['producer_plan'],plan['extraction_plan']]):assert launch['plan_sha256']==sha(path)
    groups=grouped_jobs(json.loads(Path(producer['jobs']).read_text()));assert len(groups)==405
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binding=out/'run_plan.json'
    if binding.exists():assert binding.read_bytes()==args.plan.read_bytes()
    else:binding.write_bytes(args.plan.read_bytes())
    completed={};pending={};submitted=set()
    for group in groups:
        dp=out/group/'disposition.json'
        if not dp.exists():continue
        r=json.loads(dp.read_text());assert r['group']==group
        ap=Path(r['attempt_receipt']);assert sha(ap)==r['attempt_receipt_sha256']
        for name,h in json.loads(ap.read_text())['artifacts'].items():assert sha(ap.parent/name)==h
        if 'receipt' in r:assert sha(r['receipt'])==r['receipt_sha256']
        completed[group]=r;submitted.add(group)
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        while True:
            verify();states=[dependency_state(launch) for launch in launches]
            if any(s is None for s in states):time.sleep(1);continue
            for future,group in list(pending.items()):
                if not future.done():continue
                completed[group]=future.result();del pending[future]
                print('Categorical quartet dispositions',len(completed),'/405',group,completed[group]['status'],flush=True)
            for group,rows in groups.items():
                if len(pending)>=plan['workers']:break
                if group in submitted or not ready(rows,extraction['output']):continue
                assert shutil.disk_usage(out).free>=plan['minimum_free_disk_bytes']
                pending[pool.submit(process_group,group,plan)]=group;submitted.add(group)
            status=dict(completed=len(completed),running=len(pending),pending_groups=len(groups)-len(submitted),dependency_states=states)
            write_json(out/'progress.json',status)
            terminal=all(s['ActiveState'] in ['inactive','failed'] for s in states)
            if terminal and not pending:break
            time.sleep(30)
    verify()
    write_json(out/'receipt.json',dict(status='terminal_full_categorical_report_accounting_requires_review',plan_sha256=ph,dependency_states=states,
        groups={group:completed.get(group,dict(status='unresolved_no_verified_extracted_quartet')) for group in groups},
        completed_groups=len(completed),unresolved_groups=len(groups)-len(completed),
        status_counts=dict(Counter(r['status'] for r in completed.values())),
        scope='All 405 groups accounted. Completed marginal categorical reports do not qualify joint posteriors or ancestral structures. Missing and failed groups remain explicit.'))


if __name__=='__main__':main()
