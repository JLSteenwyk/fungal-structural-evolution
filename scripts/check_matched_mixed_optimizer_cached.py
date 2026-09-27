#!/usr/bin/env python3
"""Check all-face optimizer against direct dense likelihood and Powell search."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from fit_matched_mixed_model_cached import fit_variance_ratios
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng=np.random.default_rng(81729)
    bg=np.repeat(np.arange(24),4)
    family=bg//3
    n=len(bg)
    z=(bg[:,None]==np.arange(24)[None,:]).astype(float)
    h=(family[:,None]==np.arange(8)[None,:]).astype(float)
    factor=rng.normal(size=(n,4))/2
    x=np.column_stack([np.ones(n),rng.normal(size=(n,2))])
    kernels=[z@z.T,h@h.T,factor@factor.T]
    out=Path('results/model_validation/matched-mixed-optimizer-cached-20260927-v1')
    out.mkdir(parents=True,exist_ok=False)
    summaries=[]
    for label,ratios in [('residual_only',[0.,0.,0.]),('mixed',[.7,1.3,.4]),('strong_background',[10.,.2,.1])]:
        covariance=np.eye(n)+sum(v*k for v,k in zip(ratios,kernels))
        y=x@np.array([.2,-.1,.3])+np.linalg.cholesky(covariance)@rng.normal(size=n)
        fast=fit_variance_ratios(bg,family,factor,x,y)
        def dense(parameters,details=False):
            v=np.eye(n)+sum(r*k for r,k in zip(np.expm1(parameters),kernels))
            sign,logdet=np.linalg.slogdet(v); assert sign==1
            inverse_x=np.linalg.solve(v,x)
            information=x.T@inverse_x
            beta=np.linalg.solve(information,x.T@np.linalg.solve(v,y))
            residual=y-x@beta
            q=residual@np.linalg.solve(v,residual)
            scale=q/(n-x.shape[1])
            objective=.5*(logdet+np.linalg.slogdet(information)[1]+(n-x.shape[1])*(1+np.log(2*np.pi*scale)))
            return (float(objective),beta,float(scale)) if details else float(objective)
        for candidate in fast['candidates']:
            np.testing.assert_allclose(candidate['objective'],dense(candidate['theta']),rtol=1e-9,atol=1e-8)
        direct,beta,scale=dense(fast['log1p_ratios'],True)
        np.testing.assert_allclose(fast['beta'],beta,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose(fast['profiled_scale'],scale,rtol=1e-8,atol=1e-9)
        references=[]
        for initial in [np.zeros(3),np.ones(3),np.array(fast['log1p_ratios'])]:
            ref=minimize(dense,initial,method='Powell',bounds=[(0,np.log1p(10000.))]*3,
                         options=dict(maxiter=500,ftol=1e-11,xtol=1e-9))
            references.append(dict(objective=float(ref.fun),theta=ref.x.tolist(),success=bool(ref.success),message=str(ref.message)))
        best_reference=min(r['objective'] for r in references)
        assert abs(direct-best_reference)<1e-5,(label,direct,best_reference)
        original_root=Path('results/model_validation/matched-mixed-optimizer-20260927-v1')
        original_receipt=json.loads((original_root/'receipt.json').read_text())
        assert sha(original_root/(label+'.json'))==original_receipt['artifacts'][label+'.json']
        original=json.loads((original_root/(label+'.json')).read_text())['fit']
        np.testing.assert_allclose(fast['negative_profiled_reml'],original['negative_profiled_reml'],rtol=1e-9,atol=1e-6)
        np.testing.assert_allclose(fast['beta'],original['beta'],rtol=1e-4,atol=1e-5)
        assert fast['zero_boundary']==original['zero_boundary']
        assert fast['status']==original['status']
        assert len(fast['candidates'])==22
        assert [(c['active'],c['start_ratio']) for c in fast['candidates']]==[(c['active'],c['start_ratio']) for c in original['candidates']]
        # Upper-bound clipping and genuine zero-component optima remain explicit.
        result=dict(label=label,generating_ratios=ratios,fit=fast,dense_powell_references=references,
                    absolute_best_objective_difference=abs(direct-best_reference))
        (out/(label+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        summaries.append(dict(label=label,status=fast['status'],ratios=fast['ratios'],
                              objective_difference=abs(direct-best_reference),attempts=len(fast['candidates']),wall_seconds=fast['wall_seconds'],evaluation_methods=fast['evaluation_methods']))
        print(json.dumps(summaries[-1]),flush=True)
    result=dict(status='passed_cached_optimizer_dense_and_original_solution_checks',cases=summaries,
                source_pins={p:sha(p) for p in ['scripts/matched_mixed_covariance.py','scripts/fit_matched_mixed_model_cached.py','scripts/cached_matched_likelihood.py','results/model_validation/matched-mixed-optimizer-20260927-v1/receipt.json',__file__]},
                artifacts={p.name:sha(p) for p in out.iterdir()},
                scope='All 22 boundary/start candidate likelihoods per synthetic case checked by dense algebra; best objective independently compared with three dense Powell starts. Cached solutions also checked against original direct optimizer objectives, coefficients, boundaries, statuses and complete start/face enumeration. Checks computational optimization only, not recovery of generating parameters, interval coverage, biological adequacy or a project-data pilot.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    Path('metadata/matched_mixed_optimizer_cached_checks_20260927.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
