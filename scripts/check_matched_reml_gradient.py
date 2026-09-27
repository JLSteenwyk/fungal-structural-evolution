#!/usr/bin/env python3
"""Verify analytic derivatives using an independent dense REML score identity."""
import itertools
import json
from pathlib import Path
import numpy as np
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng=np.random.default_rng(74921);n=60
    bg=np.repeat(np.arange(15),4);family=bg//3
    z=(bg[:,None]==np.arange(15)[None,:]).astype(float)
    h=(family[:,None]==np.arange(5)[None,:]).astype(float)
    x=np.column_stack([np.ones(n),rng.normal(size=(n,3))]);y=rng.normal(size=n)
    raw=rng.normal(size=(n,4));checked=0;maximum=0.
    for factor in [np.zeros((n,0)),raw,np.column_stack([raw,raw[:,0]])]:
        cache=CachedMatchedLikelihood(bg,family,factor,x,y)
        kernels=[z@z.T,h@h.T,factor@factor.T]
        for ratios in itertools.product([0.,.3,10000.],repeat=3):
            actual=evaluate_gradient(cache,ratios)
            v=np.eye(n)+sum(r*k for r,k in zip(ratios,kernels))
            inverse=np.linalg.solve(v,np.eye(n));inverse_x=inverse@x
            projection=inverse-inverse_x@np.linalg.solve(x.T@inverse_x,inverse_x.T)
            py=projection@y;q=float(y@py)
            expected=np.array([.5*(np.trace(projection@k)-(n-x.shape[1])*(py@k@py)/q) for k in kernels])
            np.testing.assert_allclose(actual['ratio_gradient'],expected,rtol=1e-7,atol=1e-7)
            np.testing.assert_allclose(actual['log1p_ratio_gradient'],expected*(1+np.array(ratios)),rtol=1e-7,atol=1e-7)
            np.testing.assert_allclose(actual['negative_profiled_reml'],cache.evaluate(*ratios)['negative_profiled_reml'],rtol=1e-9,atol=1e-9)
            maximum=max(maximum,float(np.max(abs(actual['ratio_gradient']-expected))));checked+=1
    result=dict(status='passed_analytic_matched_reml_gradient_dense_score_checks',cases=checked,ratio_derivatives_checked=checked*3,
                maximum_absolute_ratio_gradient_error=maximum,pins={p:sha(p) for p in ['scripts/matched_reml_gradient.py','scripts/cached_matched_likelihood.py',__file__]},
                scope='Independent dense residual-projector trace/score identity checks every component derivative across zero/moderate/large ratios and empty/duplicated species factors. Diagnostic implementation only; running optimizer and flags unchanged.')
    Path('metadata/matched_reml_gradient_checks_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
