"""Continue every diagnosed premature full-face endpoint with fresh solver state.

Diagnostic repair experiment only; original refits and coverage remain intact.
"""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from simulate_matched_working_model import simulate


def main():
    threadpool_limits(1)
    source=Path('results/model_validation/matched-kr-start-diagnostics-20260928-v1')
    receipt=json.loads((source/'receipt.json').read_text())
    for name,digest in receipt['artifacts'].items():assert sha(source/name)==digest
    frame=pd.read_csv(source/'full_face_candidates.tsv',sep='\t',dtype={'start':str})
    targets=frame[frame.diagnosis=='higher_objective_nonstationary']
    assert len(targets)==185
    out=Path('results/model_validation/matched-kr-start-retries-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    plan=dict(source_receipt_sha256=sha(source/'receipt.json'),targets=185,maximum_rounds=3,
        projected_gradient_tolerance=1e-6,method='L-BFGS-B',maxiter=1000,ftol=1e-14,gtol=1e-8,maxls=80,
        resources=dict(cpus=1,planning_minutes=[1,30],output_gib=.1,gpu=False,paid_cost=0),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_reml_gradient.py'),
             Path('scripts/cached_matched_likelihood.py'),Path('scripts/matched_mixed_covariance.py'),
             Path('scripts/simulate_matched_working_model.py')]})
    write_json(out/'run_plan.json',plan);ph=sha(out/'run_plan.json');rows=[];started=time.perf_counter()
    for index,row in enumerate(targets.itertuples(index=False)):
        assert sha(row.source_path)==row.source_sha256
        saved=json.loads(Path(row.source_path).read_text());case=saved['case'];f=saved['refit']['fit']
        with np.load(case['design'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        y=simulate(a['background'],a['family'],a['factor'],a['design'],a['beta'],case['scale'],case['ratios'],1,
            np.random.default_rng(np.random.SeedSequence(saved['refit']['seed_entropy'])))[:,0]
        assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest()==saved['refit']['response_sha256']
        cache=CachedMatchedLikelihood(a['background'],a['family'],a['factor'],a['design'],y)
        old=next(c for c in f['candidates'] if c['start']==row.start and all(c['active']))
        theta=np.array(old['theta']);upper=np.log1p(f['maximum_ratio']);attempts=[]
        def objective(t):
            v=evaluate_gradient(cache,np.expm1(t));return v['negative_profiled_reml'],v['log1p_ratio_gradient']
        for round_index in range(plan['maximum_rounds']):
            fit=minimize(objective,theta,jac=True,method=plan['method'],bounds=[(0,upper)]*3,
                options={k:plan[k] for k in ['maxiter','ftol','gtol','maxls']})
            value,g=objective(fit.x)
            pg=np.where(fit.x<=1e-7,np.minimum(g,0),np.where(fit.x>=upper-1e-7,np.maximum(g,0),g))
            direct=profiled_reml(MatchedCovariance(a['background'],a['family'],a['factor'],1.,*np.expm1(fit.x)),a['design'],y)
            np.testing.assert_allclose(value,direct['negative_profiled_reml'],rtol=1e-9,atol=1e-7)
            attempts.append(dict(round=round_index,theta=fit.x.tolist(),objective=float(value),
                optimizer_success=bool(fit.success),message=str(fit.message),evaluations=int(fit.nfev),
                projected_gradient_maximum=float(np.max(abs(pg))),
                direct_objective_error=abs(value-direct['negative_profiled_reml'])))
            assert value<=objective(theta)[0]+1e-7,'Restart worsened objective'
            theta=fit.x
            if np.max(abs(pg))<=plan['projected_gradient_tolerance']:break
        final=attempts[-1]
        result=dict(plan_sha256=ph,case=row.case,replicate=int(row.replicate),start=row.start,
            source_path=row.source_path,source_sha256=row.source_sha256,original_candidate=old,
            attempts=attempts,original_selected_objective=f['negative_profiled_reml'],
            final_minus_original_selected=final['objective']-f['negative_profiled_reml'],
            stationarity_pass=final['projected_gradient_maximum']<=plan['projected_gradient_tolerance'])
        path=out/f'retry-{index:04d}.json';write_json(path,result)
        rows.append(dict(path=path.name,sha256=sha(path),stationarity_pass=result['stationarity_pass'],
                         final_minus_original_selected=result['final_minus_original_selected'],rounds=len(attempts)))
    result=dict(status='completed_all_185_nonstationary_start_retry_dispositions',plan_sha256=ph,
        source_receipt_sha256=sha(source/'receipt.json'),retries=len(rows),
        stationary=sum(r['stationarity_pass'] for r in rows),
        final_agrees_with_original_selected=sum(abs(r['final_minus_original_selected'])<=1e-5 for r in rows),
        final_improves_original_selected=sum(r['final_minus_original_selected'] < -1e-5 for r in rows),
        elapsed_seconds=time.perf_counter()-started,rows=rows,
        scope='Endpoint continuation experiment; all original results preserved. No automatic '
              'replacement, coverage change, repaired-fit qualification or global optimality claim.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/matched_kr_start_retries_20260928.json'),result)
    print({k:v for k,v in result.items() if k!='rows'})


if __name__=='__main__':main()
