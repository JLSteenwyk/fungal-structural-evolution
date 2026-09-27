#!/usr/bin/env python3
"""Verify cached likelihood against checked direct residual calculations."""
import itertools
import json
import time
from pathlib import Path
import numpy as np
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng=np.random.default_rng(83157)
    bg=np.repeat(np.arange(18),5);family=bg//3;n=len(bg)
    x=np.column_stack([np.ones(n),rng.normal(size=(n,3))]);y=rng.normal(size=n)
    raw=rng.normal(size=(n,5))
    cases=0;maximum=0.;max_beta=0.
    for factor in [np.zeros((n,0)),raw,np.column_stack([raw,raw[:,0],np.zeros(n)])]:
        cache=CachedMatchedLikelihood(bg,family,factor,x,y)
        for b,f,p in itertools.product([0.,.2,10000.],[0.,.5,10000.],[0.,.3,10000.]):
            actual=cache.evaluate(b,f,p)
            expected=profiled_reml(MatchedCovariance(bg,family,factor,1.,b,f,p),x,y)
            extreme=max(b,f,p)>=10000
            for key in ['negative_profiled_reml','beta','profiled_scale','residual_quadratic','conditional_beta_covariance']:
                np.testing.assert_allclose(actual[key],expected[key],rtol=1e-7 if extreme else 1e-9,atol=1e-7 if extreme else 1e-10)
            maximum=max(maximum,abs(actual['negative_profiled_reml']-expected['negative_profiled_reml']))
            max_beta=max(max_beta,float(np.max(abs(actual['beta']-expected['beta']))))
            cases+=1
    # Deliberately cancellation-prone response; fallback must preserve the checked path.
    large=x@np.array([1e7,-3e6,2e6,4e6])+rng.normal(size=n)
    cache=CachedMatchedLikelihood(bg,family,raw,x,large)
    actual=cache.evaluate(.2,.5,.3)
    expected=profiled_reml(MatchedCovariance(bg,family,raw,1.,.2,.5,.3),x,large)
    assert actual['evaluation_method']=='direct_residual_fallback'
    for key in ['negative_profiled_reml','beta','profiled_scale','conditional_beta_covariance']:
        np.testing.assert_array_equal(actual[key],expected[key])
    # Full-size numerical comparison and timing, with unchanged parameters/results.
    timing=[]
    for n,g,nf in [(148,120,44),(2320,1500,400),(10963,6478,1634)]:
        bg=np.concatenate([np.arange(g),rng.integers(0,g,n-g)])
        gf=np.concatenate([np.arange(nf),rng.integers(0,nf,g-nf)])
        family=gf[bg];factor=rng.normal(size=(n,242))/np.sqrt(242)
        x=np.column_stack([np.ones(n),rng.normal(size=(n,4))]);y=rng.normal(size=n)
        before=time.perf_counter();cache=CachedMatchedLikelihood(bg,family,factor,x,y)
        preparation=time.perf_counter()-before
        fast_times=[];direct_times=[]
        for repeat in range(5):
            before=time.perf_counter();actual=cache.evaluate(.2,.5,.3);fast_times.append(time.perf_counter()-before)
            before=time.perf_counter();expected=profiled_reml(MatchedCovariance(bg,family,factor,1.,.2,.5,.3),x,y);direct_times.append(time.perf_counter()-before)
            for key in ['negative_profiled_reml','beta','profiled_scale','conditional_beta_covariance']:
                np.testing.assert_allclose(actual[key],expected[key],rtol=1e-9,atol=1e-9)
        timing.append(dict(records=n,backgrounds=g,families=nf,preparation_seconds=preparation,
                           cached_seconds=fast_times,direct_seconds=direct_times,
                           median_speedup=float(np.median(direct_times)/np.median(fast_times))))
        print(json.dumps(timing[-1]),flush=True)
    result=dict(status='passed_cached_matched_likelihood_reference_checks',boundary_factor_cases=cases,
                cancellation_fallback_cases=1,maximum_objective_error=maximum,maximum_coefficient_error=max_beta,
                full_size_timing=timing,pins={p:sha(p) for p in ['scripts/cached_matched_likelihood.py','scripts/matched_mixed_covariance.py',__file__]},
                scope='Full zero/high-ratio component grid and empty/duplicated species factors compared with checked direct evaluator; cancellation fixture must use direct residual fallback. Full-size synthetic timings include five paired repeats and separate setup. No model fit, optimizer equivalence, biological subsampling or full-grid ETA.')
    Path('metadata/cached_matched_likelihood_checks_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['full_size_timing','pins']},indent=2))


if __name__=='__main__':main()
