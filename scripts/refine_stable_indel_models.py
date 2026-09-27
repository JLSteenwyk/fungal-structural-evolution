#!/usr/bin/env python3
"""Refine every retained indel solution with direct event ascertainment."""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from stable_indel_likelihood import StableLikelihood
from prepare_case_ancestral_neighborhoods import sha


def fit(job):
    with threadpool_limits(limits=1):
        return run(job)


def run(job):
    for p,h in job['pins'].items():
        assert sha(p)==h,p
    source=json.loads(Path(job['source_receipt']).read_text())
    out=Path(job['output']);out.mkdir(parents=True,exist_ok=False)
    result=dict(job=job)
    if not job['character_count']:
        result['status']='no_coded_characters_no_inference'
    else:
        seq={r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
        model=StableLikelihood(Phylo.read(job['tree'],'newick'),seq,job['correction'])
        bounds=[(-12,12),(-12,12),(math.log(.005),math.log(100))]
        records=[]
        def objective(v):
            value=model.evaluate(v)
            return -value if math.isfinite(value) else 1e100
        def candidate(x,stage,optimizer=None):
            return dict(stage=stage,log_parameters=np.asarray(x).tolist(),
                parameters=np.exp(x).tolist(),log_likelihood=model.evaluate(x),
                optimizer_success=bool(optimizer.success) if optimizer is not None else None,
                optimizer_message=str(optimizer.message) if optimizer is not None else 'retained starting point',
                evaluations=int(optimizer.nfev) if optimizer is not None else 1)
        for index,old in enumerate(source['starts']):
            initial=candidate(np.array(old['log_parameters']),'initial_stable_replay')
            powell=minimize(objective,initial['log_parameters'],method='Powell',bounds=bounds,
                options=dict(ftol=1e-10,xtol=1e-6,maxiter=300))
            trial=candidate(powell.x,'Powell',powell)
            seed=max([initial,trial],key=lambda r:r['log_likelihood'])
            polish=minimize(objective,seed['log_parameters'],method='L-BFGS-B',jac='3-point',bounds=bounds,
                options=dict(ftol=1e-12,gtol=1e-6,maxiter=500,maxls=40,finite_diff_rel_step=1e-5))
            polished=candidate(polish.x,'L-BFGS-B_3point',polish)
            candidates=[initial,trial,polished]
            best=max(candidates,key=lambda r:r['log_likelihood'])
            x=np.array(best['log_parameters']);directions=[]
            for step in [1e-4,1e-3]:
                for axis in range(3):
                    for sign in [-1,1]:
                        v=x.copy();v[axis]=np.clip(v[axis]+sign*step,*bounds[axis])
                        if v[axis]==x[axis]:
                            continue
                        directions.append(dict(axis=axis,requested_step=sign*step,
                            actual_step=float(v[axis]-x[axis]),
                            likelihood_change=model.evaluate(v)-best['log_likelihood']))
            records.append(dict(source_start=index,source_reported_log_likelihood=old['log_likelihood'],
                initial_stable_log_likelihood=initial['log_likelihood'],candidates=candidates,best=best,
                signed_improvement=best['log_likelihood']-initial['log_likelihood'],
                near_bound=[bool(min(abs(v-lo),abs(v-hi))<1e-4) for v,(lo,hi) in zip(x,bounds)],
                feasible_coordinate_checks=directions))
        best=max(records,key=lambda r:r['best']['log_likelihood'])
        result.update(status='stable_refinements_produced_pending_independent_audit',starts=records,
            best=best['best'],selected_source_start=best['source_start'],log_parameter_bounds=bounds,
            refined_start_likelihood_spread=max(r['best']['log_likelihood'] for r in records)-min(r['best']['log_likelihood'] for r in records))
    (out/'receipt.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return dict(job_id=job['job_id'],status=result['status'],receipt_sha256=sha(out/'receipt.json'))


def main():
    pp=Path('metadata/stable_indel_refinement_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    results=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for result in pool.map(fit,plan['jobs']):
            results.append(result);print(result['job_id'],result['status'],flush=True)
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    (out/'receipt.json').write_text(json.dumps(dict(status='all_stable_refinements_produced_pending_audit',
        jobs=results,plan_sha256=sha(pp),scope=plan['scope']),indent=2)+'\n')


if __name__=='__main__':
    main()
