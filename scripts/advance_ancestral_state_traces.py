#!/usr/bin/env python3
"""Advance state extraction across the full immutable production chain grid."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import psutil
from ancestral_chain_attempt import run_attempt,sha,write_json


def snapshot_row(job,audit_path,iterations):
    chain,config=job['chain'],job['config']
    audit=json.loads(audit_path.read_text())
    assert audit['status']=='all_saved_alignments_and_candidate_nodes_checked' and audit['iterations']==iterations
    receipt_path=Path(audit['attempt_receipt']);assert sha(receipt_path)==audit['attempt_receipt_sha256']
    receipt=json.loads(receipt_path.read_text());assert receipt['exit_code']==0
    assert json.loads((audit_path.parent/'configuration.json').read_text())==config
    digest=hashlib.sha256(json.dumps(config,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    assert receipt['configuration_sha256']==digest
    for p,h in config['pins'].items():assert sha(p)==h
    for name,h in receipt['artifacts'].items():assert sha(receipt_path.parent/name)==h
    assert sha(audit['scalar_log'])==audit['scalar_log_sha256']
    assert chain['seed']==config['seed']
    identity=json.loads((receipt_path.parent/'process.json').read_text())
    try:
        p=psutil.Process(identity['pid'])
        assert p.create_time()!=identity['created'] or p.status()==psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:pass
    assert len(audit['candidate_samples'])==4*(iterations//10+1)
    return dict(chain=chain['chain_id'],seed=chain['seed'],audit=str(audit_path.resolve()),
                audit_sha256=sha(audit_path),saved_alignments=iterations//10+1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());digest=sha(args.plan)
    def verify():
        assert sha(args.plan)==digest
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify()
    producer=json.loads(Path(plan['producer_plan']).read_text());launch=json.loads(Path(plan['producer_launch']).read_text())
    assert launch['plan_sha256']==sha(plan['producer_plan'])
    jobs=json.loads(Path(producer['jobs']).read_text());assert len(jobs)==len({j['chain']['chain_id'] for j in jobs})==1620
    out=Path(plan['output']).resolve();out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binding=out/'run_plan.json'
    if binding.exists():assert binding.read_bytes()==args.plan.read_bytes()
    else:binding.write_bytes(args.plan.read_bytes())
    completed={}
    while True:
        verify()
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        terminal=state['ActiveState'] in ['inactive','failed']
        if not terminal:
            assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
            try:
                p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
            except psutil.NoSuchProcess:
                time.sleep(1);continue
        for job in jobs:
            chain_id=job['chain']['chain_id']
            if chain_id in completed:continue
            source=Path(producer['output'])/chain_id
            audits=sorted(source.glob('attempt-*-sample-audit.json'))
            if not audits:continue
            assert len(audits)==1, 'Multiple audited attempts require explicit review'
            row=snapshot_row(job,audits[0],producer['iterations'])
            assert json.loads(audits[0].read_text())['mapping_sha256']==sha(producer['mapping'])
            target=out/chain_id;target.mkdir(exist_ok=True)
            snapshot=target/'snapshot.json'
            payload=dict(status='verified_frozen_terminal_chain_subset',rows=[row])
            if snapshot.exists():assert json.loads(snapshot.read_text())==payload
            else:write_json(snapshot,payload)
            assert shutil.disk_usage(out).free>=plan['minimum_free_disk_bytes']
            config=dict(command=[str(Path(plan['python']).absolute()),str(Path(plan['script']).resolve()),
                '--snapshot',str(snapshot),'--inputs',str(Path(plan['inputs']).resolve()),'--output','{attempt}/states'],
                timeout_seconds=plan['per_chain_timeout_seconds'],pins={**{str(Path(p).absolute()):h for p,h in plan['pins'].items()},str(snapshot):sha(snapshot)},
                chain_id=chain_id,source_audit_sha256=sha(audits[0]))
            # run_attempt requires absolute executable paths in its pin map.
            config['pins'][config['command'][0]]=sha(config['command'][0])
            attempt=run_attempt(target/'processing',config);attempt_receipt=json.loads(attempt.read_text())
            if attempt_receipt['exit_code']!=0:
                result=dict(status='state_extraction_attempt_failed_requires_review',receipt=str(attempt),receipt_sha256=sha(attempt))
            else:
                rp=attempt.parent/'states/receipt.json';report=json.loads(rp.read_text())
                assert report['status']=='complete_independently_checked_anchored_state_traces' and report['chains']==1
                assert report['summaries'][0]['chain']==chain_id
                for summary in report['summaries']:
                    for p,h in summary['artifacts'].items():assert sha(p)==h
                for p,h in report['pins'].items():assert sha(p)==h
                result=dict(status='state_trace_complete_not_posterior_qualification',receipt=str(rp),receipt_sha256=sha(rp),attempt_receipt=str(attempt),attempt_receipt_sha256=sha(attempt),state_observations=report['state_observations'])
            write_json(target/'disposition.json',result);completed[chain_id]=result
            print('State-trace dispositions',len(completed),'/1620',chain_id,result['status'],flush=True)
        if terminal:break
        time.sleep(30)
    verify()
    write_json(out/'receipt.json',dict(status='terminal_full_state_trace_accounting_requires_review',plan_sha256=digest,producer_state=state,
        chains={j['chain']['chain_id']:completed.get(j['chain']['chain_id'],dict(status='unresolved_no_verified_terminal_chain')) for j in jobs},
        processed_chains=len(completed),unresolved_chains=1620-len(completed),
        scope='Full grid preserved with explicit unresolved/failed chains. State arrays are diagnostic inputs, not qualified posteriors.'))


if __name__=='__main__':main()
