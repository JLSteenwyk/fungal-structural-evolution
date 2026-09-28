#!/usr/bin/env python3
"""Independent dense likelihood and derivative-free checks of analytic refinement."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from refine_matched_reml_analytic import refine
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    threadpool_limits(1)
    rng = np.random.default_rng(81729)
    bg = np.repeat(np.arange(24), 4)
    family = bg//3
    n = len(bg)
    z = (bg[:,None] == np.arange(24)).astype(float)
    h = (family[:,None] == np.arange(8)).astype(float)
    factor = rng.normal(size=(n,4))/2
    x = np.column_stack([np.ones(n), rng.normal(size=(n,2))])
    summaries = []
    for label, ratios in [('residual_only',[0.,0.,0.]), ('mixed',[.7,1.3,.4]), ('strong_background',[10.,.2,.1])]:
        kernels = [z@z.T, h@h.T, factor@factor.T]
        v = np.eye(n)+sum(r*k for r,k in zip(ratios,kernels))
        y = x@np.array([.2,-.1,.3])+np.linalg.cholesky(v)@rng.normal(size=n)
        old_root = Path('results/model_validation/matched-mixed-optimizer-cached-20260927-v1')
        old_receipt = json.loads((old_root/'receipt.json').read_text())
        path = old_root/(label+'.json')
        assert sha(path) == old_receipt['artifacts'][path.name]
        old = json.loads(path.read_text())['fit']
        result = refine(bg, family, factor, x, y, old['log1p_ratios'])

        def dense(theta, details=False):
            covariance = np.eye(n)+sum(r*k for r,k in zip(np.expm1(theta),kernels))
            inv_x = np.linalg.solve(covariance,x)
            info = x.T@inv_x
            beta = np.linalg.solve(info,x.T@np.linalg.solve(covariance,y))
            residual = y-x@beta
            scale = (residual@np.linalg.solve(covariance,residual))/(n-x.shape[1])
            value = .5*(np.linalg.slogdet(covariance)[1]+np.linalg.slogdet(info)[1]+(n-x.shape[1])*(1+np.log(2*np.pi*scale)))
            return (float(value),beta,float(scale)) if details else float(value)

        assert len(result['candidates']) == 24
        assert len({tuple(c['active']) for c in result['candidates']}) == 8
        for candidate in result['candidates']:
            np.testing.assert_allclose(candidate['objective'],dense(candidate['theta']),rtol=1e-9,atol=1e-8)
        objective,beta,scale = dense(result['log1p_ratios'],True)
        np.testing.assert_allclose(result['beta'],beta,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose(result['profiled_scale'],scale,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose(result['original_objective'],old['negative_profiled_reml'],rtol=1e-9,atol=1e-8)
        references = []
        for start in [np.zeros(3), np.ones(3), np.array(result['log1p_ratios'])]:
            fit = minimize(dense,start,method='Powell',bounds=[(0.,np.log1p(10000.))]*3,
                           options=dict(maxiter=500,ftol=1e-11,xtol=1e-9))
            references.append(dict(objective=float(fit.fun),success=bool(fit.success)))
        difference = abs(objective-min(r['objective'] for r in references))
        assert difference < 1e-5, (label,difference)
        assert result['checks']['original_objective_not_worsened']
        assert result['checks']['projected_gradient_pass']
        (args.output/(label+'.json')).write_text(json.dumps(dict(fit=result,dense_powell=references,source_sha256=sha(path)),indent=2)+'\n')
        summaries.append(dict(label=label,status=result['status'],wall_seconds=result['wall_seconds'],
                              objective_difference=difference,improvement=result['objective_improvement']))
    for bad in [[-1,0,0],[0,0,float('nan')],[0,0,100],[0,0]]:
        try:
            refine(bg,family,factor,x,y,bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid initial parameters accepted')
    pins = ['scripts/refine_matched_reml_analytic.py',__file__,'scripts/matched_reml_gradient.py',
            'scripts/cached_matched_likelihood.py','scripts/matched_mixed_covariance.py',str(old_root/'receipt.json')]
    receipt = dict(status='passed_dense_and_derivative_free_refinement_checks',cases=summaries,
                   candidate_likelihoods_checked=72,invalid_inputs_rejected=4,
                   pins={p:sha(p) for p in pins},artifacts={p.name:sha(p) for p in args.output.iterdir()},
                   scope='Synthetic numerical fixtures; production review fits have not been refined by this check.')
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
