#!/usr/bin/env python3
"""Compare nested/low-rank likelihood algebra with direct dense covariance."""
import itertools
import json
from pathlib import Path
import numpy as np
from matched_mixed_covariance import MatchedCovariance, profiled_reml
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng = np.random.default_rng(20260927)
    n = 73
    bg = rng.integers(0,15,n)
    fam = bg//3
    z = (bg[:,None] == np.unique(bg)[None,:]).astype(float)
    h = (fam[:,None] == np.unique(fam)[None,:]).astype(float)
    raw = rng.normal(size=(n,5))
    factors = [np.zeros((n,0)), raw, np.column_stack([raw,raw[:,0],np.zeros(n)])]
    x = np.column_stack([np.ones(n),rng.normal(size=(n,3))])
    y = rng.normal(size=n)
    rhs = rng.normal(size=(n,9))
    worst_solve = worst_logdet = worst_likelihood = 0.
    cases = 0
    for f in factors:
        for residual,bvar,fvar,pvar in itertools.product([.3,1.,4.], [0.,.2,10.], [0.,.5,20.], [0.,.1,3.]):
            evaluator = MatchedCovariance(bg,fam,f,residual,bvar,fvar,pvar)
            dense = residual*np.eye(n)+bvar*(z@z.T)+fvar*(h@h.T)+pvar*(f@f.T)
            expected_solve = np.linalg.solve(dense,rhs)
            actual_solve = evaluator.solve(rhs)
            np.testing.assert_allclose(actual_solve,expected_solve,rtol=1e-8,atol=1e-10)
            np.testing.assert_allclose(evaluator.solve(rhs[:,0]),expected_solve[:,0],rtol=1e-8,atol=1e-10)
            sign,logdet = np.linalg.slogdet(dense)
            assert sign == 1
            np.testing.assert_allclose(evaluator.logdet,logdet,rtol=1e-10,atol=1e-10)
            inv_x = np.linalg.solve(dense,x)
            info = x.T@inv_x
            beta = np.linalg.solve(info, x.T@np.linalg.solve(dense,y))
            r = y-x@beta
            q = r@np.linalg.solve(dense,r)
            scale = q/(n-x.shape[1])
            likelihood = .5*(logdet+np.linalg.slogdet(info)[1]+(n-x.shape[1])*(1+np.log(2*np.pi*scale)))
            result = profiled_reml(evaluator,x,y)
            np.testing.assert_allclose(result['beta'],beta,rtol=1e-8,atol=1e-10)
            np.testing.assert_allclose(result['conditional_beta_covariance'],scale*np.linalg.inv(info),rtol=1e-8,atol=1e-10)
            np.testing.assert_allclose(result['profiled_scale'],scale,rtol=1e-9,atol=1e-10)
            np.testing.assert_allclose(result['negative_profiled_reml'],likelihood,rtol=1e-9,atol=1e-9)
            worst_solve=max(worst_solve,float(np.max(abs(actual_solve-expected_solve))))
            worst_logdet=max(worst_logdet,float(abs(evaluator.logdet-logdet)))
            worst_likelihood=max(worst_likelihood,float(abs(result['negative_profiled_reml']-likelihood)))
            # Row order cannot affect determinant, coefficients, or restricted likelihood.
            order=rng.permutation(n)
            perm=MatchedCovariance(bg[order],fam[order],f[order],residual,bvar,fvar,pvar)
            reordered=profiled_reml(perm,x[order],y[order])
            np.testing.assert_allclose(reordered['beta'],result['beta'],rtol=1e-8,atol=1e-10)
            np.testing.assert_allclose(reordered['negative_profiled_reml'],result['negative_profiled_reml'],rtol=1e-9,atol=1e-9)
            cases+=1
    failures=0
    bad_family=fam.copy()
    repeated=np.flatnonzero(bg==bg[0])
    assert len(repeated)>1
    bad_family[repeated[0]]=100
    invalid=[lambda:MatchedCovariance(bg,bad_family,raw),
             lambda:MatchedCovariance(bg,fam,raw,residual=0),
             lambda:MatchedCovariance(bg,fam,raw,background_variance=-1),
             lambda:MatchedCovariance(bg,fam,raw,species_variance=np.nan),
             lambda:profiled_reml(MatchedCovariance(bg,fam,raw),np.column_stack([x,x[:,0]]),y)]
    for call in invalid:
        try:
            call()
        except ValueError:
            failures+=1
        else:
            raise AssertionError('Invalid input unexpectedly accepted')
    receipt=dict(status='passed_dense_reference_matched_mixed_covariance_checks',
                 variance_factor_cases=cases,permuted_cases=cases,expected_rejections=failures,
                 maximum_absolute_solve_error=worst_solve,maximum_absolute_logdet_error=worst_logdet,
                 maximum_absolute_profiled_reml_error=worst_likelihood,
                 pins={p:sha(p) for p in ['scripts/matched_mixed_covariance.py',__file__]},
                 scope='Synthetic numerical validation across zero/nonzero components, three residual scales, empty/full/duplicated species factors and row permutations. Direct dense solves, determinants, generalized least squares and restricted likelihood agree. Not a dataset pilot, optimizer verification, confidence-interval calibration or biological result.')
    Path('metadata/matched_mixed_covariance_checks_20260927.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
