#!/usr/bin/env python3
"""Restartable complete output audit; unresolved fitting errors remain explicit."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from threadpoolctl import threadpool_limits
from check_whole_protein_fit_payload import check_payload
from screen_duplication_alignment_reuse import sha

FACTORS = {}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def initialize(root):
    global FACTORS
    threadpool_limits(limits=1)
    FACTORS = {}
    for path in sorted(Path(root).glob('*.npz')):
        with np.load(path,allow_pickle=False) as a:
            FACTORS[path.stem] = a['factor']
    assert len(FACTORS) == 5


def audit_case(task):
    item, rows, inputs, output, fit_plan_hash, audit_plan_hash = task
    identifier = item['fit_input_id']
    assert {r['tree'] for r in rows} == set(FACTORS) and len(rows) == 5
    path = Path(inputs)/item['path']
    assert sha(path) == item['sha256']
    hashes = {r['path']:r['sha256'] for r in rows}
    for p,h in hashes.items():
        assert sha(p) == h,p
    identity = dict(fit_input_id=identifier,input_sha256=item['sha256'],fit_hashes=hashes,fit_plan_sha256=fit_plan_hash,audit_plan_sha256=audit_plan_hash)
    proof_path = Path(output)/'inputs'/identifier[:2]/(identifier+'.json')
    if proof_path.exists():
        proof = json.loads(proof_path.read_text())
        assert proof['identity'] == identity and proof['results_sha256'] == digest(proof['results'])
        return proof['results']
    with np.load(path,allow_pickle=False) as a:
        matrix,bg,fam,pattern_rows,identities = a['matrix'],a['background'],a['family'],a['pattern_rows'],a['row_identity']
    spec = item['recipe']['specification']
    assert digest(spec) == identifier
    assert hashlib.sha256(matrix.tobytes()).hexdigest() == spec['values_sha256']
    assert hashlib.sha256(identities.tobytes()).hexdigest() == spec['ordered_identity_sha256']
    results = []
    for row in sorted(rows,key=lambda r:r['tree']):
        saved = json.loads(Path(row['path']).read_text())
        assert saved['plan_sha256'] == fit_plan_hash and saved['input_sha256'] == item['sha256']
        assert saved['fit_input_id'] == row['fit_input_id'] == identifier and saved['tree'] == row['tree']
        assert saved['specification'] == spec and saved['payload_sha256'] == digest(saved['payload'])
        payload = saved['payload']
        assert row['status'] == payload['status']
        if payload['status'] == 'fit_error_requires_review':
            assert payload['error_type'] and payload['traceback']
            result = dict(status=payload['status'],candidates_checked=0,numerical_fit_verified=False)
        else:
            full_factor = FACTORS[row['tree']]
            assert np.all(pattern_rows>=0) and np.all(pattern_rows<len(full_factor))
            result = check_payload(payload,bg,fam,full_factor[pattern_rows],matrix,spec['columns'])
            result['numerical_fit_verified'] = True
        results.append(dict(fit_input_id=identifier,tree=row['tree'],**result))
    proof = dict(identity=identity,results=results,results_sha256=digest(results))
    proof_path.parent.mkdir(parents=True,exist_ok=True)
    temporary = proof_path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(proof,indent=2)+'\n')
    temporary.replace(proof_path)
    return results


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    launch=json.loads(Path(plan['fit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p=psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    producer=json.loads(Path(launch['plan']).read_text());fit_ph=sha(launch['plan'])
    assert fit_ph==launch['plan_sha256']
    root=Path(producer['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_whole_protein_ml_dispositions_pending_full_audit' and receipt['plan_sha256']==fit_ph
    bindings={str(args.plan):ph,**plan['pins'],str(root/'receipt.json'):sha(root/'receipt.json')}
    bindings.update(json.loads((root/'source_bindings.json').read_text()))
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify()
    rows=defaultdict(list);seen=set()
    for line in (root/'fit_manifest.jsonl').open():
        row=json.loads(line);key=row['fit_input_id'],row['tree'];assert key not in seen;seen.add(key);rows[key[0]].append(row)
    assert len(seen)==receipt['tree_fit_dispositions']==375350 and len(rows)==receipt['unique_inputs']==75070
    output=Path(plan['output']);output.mkdir(parents=True,exist_ok=True)
    lock=(output/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    saved_plan=output/'run_plan.json'
    if saved_plan.exists():assert saved_plan.read_bytes()==args.plan.read_bytes()
    else:saved_plan.write_bytes(args.plan.read_bytes())
    pending=set();submitted=set();counts=Counter();total=candidates=0
    with (output/'audit_manifest.jsonl').open('w') as manifest,ProcessPoolExecutor(max_workers=plan['workers'],mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(producer['factors'],)) as pool:
        def collect():
            nonlocal pending,total,candidates
            done,pending=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                for result in future.result():
                    counts[result['status']]+=1;total+=1;candidates+=result['candidates_checked'];manifest.write(json.dumps(result)+'\n')
                manifest.flush()
            print('Audited whole-protein dispositions',total,'/375350',dict(counts),flush=True)
        for line in (Path(producer['inputs'])/'input_manifest.jsonl').open():
            item=json.loads(line);key=item['fit_input_id'];assert key not in submitted;submitted.add(key)
            pending.add(pool.submit(audit_case,(item,rows[key],producer['inputs'],str(output),fit_ph,ph)))
            if len(pending)>=2*plan['workers']:collect()
        while pending:collect()
    assert submitted==set(rows) and total==375350 and dict(counts)==receipt['status_counts']
    verify()
    result=dict(status='complete_full_whole_protein_ml_output_audit',tree_fit_dispositions=total,candidate_likelihoods_checked=candidates,status_counts=dict(counts),source_receipt_sha256=sha(root/'receipt.json'),audit_plan_sha256=ph,producer_terminal_state=state,artifacts={p.name:sha(p) for p in [output/'run_plan.json',output/'audit_manifest.jsonl']},scope='All dispositions accounted for. Numeric payloads replayed; error records retained with numerical_fit_verified=false. Review flags remain unresolved. No model preference, calibrated uncertainty, global optimum or biological inference.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
