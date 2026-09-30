#!/usr/bin/env python3
"""Restartable numerical follow-up of every optimization flag in the full grid."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import json
import multiprocessing
from pathlib import Path
import time
import traceback
from threadpoolctl import threadpool_limits
from whole_protein_flag_followup import sources, arrays, followup, FLAG, ERROR
from screen_duplication_alignment_reuse import sha
from run_whole_protein_ml import digest


def initialize():threadpool_limits(limits=1)


def work(task):
    item, recipe, fit_plan, output, plan_hash = task
    identity = dict(fit_input_id=item['fit_input_id'],tree=item['tree'],
        original_fit=item['path'],original_fit_sha256=item['sha256'],
        input_sha256=recipe['sha256'],plan_sha256=plan_hash)
    root = Path(output); target = root/'cases'/item['fit_input_id'][:2]/(item['fit_input_id']+'-'+item['tree']+'.json')
    if target.exists():
        saved = json.loads(target.read_text())
        assert saved['identity'] == identity and saved['result_sha256'] == digest(saved['result'])
    else:
        original, matrix, bg, family, factor, x, y, scales = arrays(item, recipe, fit_plan)
        began = time.perf_counter()
        try:
            result = followup(original,bg,family,factor,x,y,scales)
        except Exception as error:
            result = dict(status='numerical_followup_error_requires_review',
                error_type=type(error).__name__,error=str(error),traceback=traceback.format_exc())
        saved = dict(identity=identity,result=result,result_sha256=digest(result),
            original_status=original['payload']['status'],original_checks=original['payload']['checks'],
            columns=original['specification']['columns'],covariate_scales=scales.tolist(),
            records=len(matrix),wall_seconds=time.perf_counter()-began,scientific_eligibility=False)
        target.parent.mkdir(parents=True,exist_ok=True)
        temporary = target.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(saved,indent=2,allow_nan=False)+'\n');temporary.replace(target)
    return dict(fit_input_id=item['fit_input_id'],tree=item['tree'],path=str(target.relative_to(root)),
        sha256=sha(target),status=saved['result']['status'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',required=True,type=Path)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    fit_plan, production, flags, errors, recipes, bindings = sources(plan)
    bindings[str(args.plan)] = ph
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify();root=Path(plan['output']);root.mkdir(parents=True,exist_ok=True)
    lock=(root/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    for name,value in [('run_plan.json',args.plan.read_text()),('source_bindings.json',json.dumps(bindings,indent=2,sort_keys=True)+'\n')]:
        path=root/name
        if path.exists():assert path.read_text()==value
        else:path.write_text(value)
    with (root/'scope_dispositions.jsonl').open('w') as f:
        for item in production.values():
            disposition=('optimization_flag_followup_required' if item['status']==FLAG else
                'original_fit_error_retained_requires_review' if item['status']==ERROR else 'original_numerical_pass_retained')
            f.write(json.dumps(dict(**item,followup_scope_status=disposition,scientific_eligibility=False))+'\n')
    pending=set();counts=Counter();completed=0
    with (root/'followup_manifest.jsonl').open('w') as manifest, ProcessPoolExecutor(
        max_workers=plan['workers'],mp_context=multiprocessing.get_context('spawn'),initializer=initialize) as pool:
        def collect():
            nonlocal pending,completed
            done,pending=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                row=future.result();manifest.write(json.dumps(row)+'\n');counts[row['status']]+=1;completed+=1
            manifest.flush();print('Completed full-grid optimization followups',completed,'/',len(flags),dict(counts),flush=True)
        for item in flags:
            pending.add(pool.submit(work,(item,recipes[item['fit_input_id']],fit_plan,str(root),ph)))
            if len(pending)>=2*plan['workers']:collect()
        while pending:collect()
    assert completed==len(flags)
    verify()
    receipt=dict(status='complete_full_grid_optimization_followup_pending_independent_readback',
        plan_sha256=ph,full_dispositions=len(production),flagged_fits=len(flags),
        original_fit_errors_retained=len(errors),followup_status_counts=dict(counts),
        artifacts={name:sha(root/name) for name in ['run_plan.json','source_bindings.json','scope_dispositions.jsonl','followup_manifest.jsonl']},
        scientific_eligibility=False,scope='Entire audited 375350-fit grid accounted for, all optimization flags attempted with separate stricter 24-candidate all-face refinement and, when needed, five bounded full-face recoveries. Original source fits/flags unchanged. Errors and unresolved numerical review retained. Independent readback, complete comparison integration, calibration and biological adequacy still required.')
    (root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
