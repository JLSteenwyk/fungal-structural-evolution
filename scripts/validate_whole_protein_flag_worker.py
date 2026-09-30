#!/usr/bin/env python3
"""Exercise producer, saved-output checker and restart on five immutable cases."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
from run_whole_protein_flag_followup import work, initialize
from readback_whole_protein_flag_followup import work as readback
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',required=True,type=Path)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    bindings={str(args.plan):ph,**plan['pins']}
    registry=json.loads(Path(plan['registry']).read_text())
    assert registry['status']=='completed_frozen_five_fit_numerical_recovery_candidates'
    for p,h in registry['source_hashes'].items():
        assert p not in bindings or bindings[p]==h;bindings[p]=h
    def verify():
        for p,h in bindings.items():assert sha(p)==h,p
    verify();fit_plan=json.loads(Path(plan['fit_plan']).read_text())
    wanted={r['fit_input_id'] for r in registry['candidates']};recipes={}
    with (Path(fit_plan['inputs'])/'input_manifest.jsonl').open() as f:
        for line in f:
            row=json.loads(line)
            if row['fit_input_id'] in wanted:recipes[row['fit_input_id']]=row
    assert len(recipes)==len(wanted)==5
    root=Path(plan['output']);root.mkdir(parents=True,exist_ok=False)
    tasks=[]
    for entry in registry['candidates']:
        item=dict(fit_input_id=entry['fit_input_id'],tree=entry['tree'],
            path=entry['original_production_fit'],sha256=entry['original_production_sha256'])
        tasks.append((item,recipes[item['fit_input_id']],fit_plan,str(root),ph))
    with ProcessPoolExecutor(max_workers=plan['workers'],mp_context=multiprocessing.get_context('spawn'),initializer=initialize) as pool:
        rows=list(pool.map(work,tasks,chunksize=1))
    assert len(rows)==5
    fit_plan['followup_plan_sha256']=ph;checked=[]
    for row,task in zip(rows,tasks):
        item,recipe,_,_,_=task
        result=readback((row,item,recipe,fit_plan,str(root)))
        assert result['numerical_result_verified']
        assert result['status']=='numerical_followup_passed_pending_independent_readback'
        before=(root/row['path']).read_bytes();resumed=work(task)
        assert resumed==row and (root/row['path']).read_bytes()==before
        checked.append(result);print('Producer/readback/restart validated',len(checked),'/5',flush=True)
    verify()
    proof=dict(status='passed_actual_five_case_flag_worker_readback_and_restart_validation',
        cases=5,results=checked,byte_identical_restarts=True,source_hashes=bindings,
        dispositions=rows,scope='All five previously audited immutable source cases exercised through the actual producer worker, separate serialized checker and identical-plan restart. This is implementation verification within the full workflow, not a pilot or full-grid numerical/scientific acceptance.')
    with (root/'receipt.json').open('x') as f:json.dump(proof,f,indent=2);f.write('\n')


if __name__=='__main__':main()
