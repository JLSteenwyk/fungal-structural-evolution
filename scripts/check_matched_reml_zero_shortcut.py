"""Verify exact zero-component shortcut with dense scores and original evaluator."""
import itertools
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient as original
from matched_reml_gradient_fast import evaluate_gradient as fast


def main():
    threadpool_limits(1);rng=np.random.default_rng(20260928);n=60
    bg=np.repeat(np.arange(15),4);fam=bg//3
    x=np.column_stack([np.ones(n),rng.normal(size=(n,3))]);y=rng.normal(size=n)
    raw=rng.normal(size=(n,4));cases=0;maximum=0.
    for f in [np.zeros((n,0)),raw,np.column_stack([raw,raw[:,0]])]:
        cache=CachedMatchedLikelihood(bg,fam,f,x,y)
        kernels=[(bg[:,None]==bg[None,:]).astype(float),(fam[:,None]==fam[None,:]).astype(float),f@f.T]
        for ratios in itertools.product([0.,.3,10000.],repeat=3):
            result=fast(cache,ratios);old=original(cache,ratios)
            np.testing.assert_allclose(result['ratio_gradient'],old['ratio_gradient'],rtol=1e-12,atol=1e-10)
            np.testing.assert_allclose(result['negative_profiled_reml'],old['negative_profiled_reml'],rtol=1e-12,atol=1e-10)
            v=np.eye(n)+sum(r*k for r,k in zip(ratios,kernels));inverse=np.linalg.solve(v,np.eye(n));ix=inverse@x
            projection=inverse-ix@np.linalg.solve(x.T@ix,ix.T);py=projection@y;q=y@py
            score=np.array([.5*(np.trace(projection@k)-(n-x.shape[1])*(py@k@py)/q) for k in kernels])
            np.testing.assert_allclose(result['ratio_gradient'],score,rtol=1e-7,atol=1e-7)
            maximum=max(maximum,float(np.max(abs(result['ratio_gradient']-score))));cases+=1
    n=600;bg=np.repeat(np.arange(100),6);fam=bg//5
    cache=CachedMatchedLikelihood(bg,fam,rng.normal(size=(n,250))/16,np.column_stack([np.ones(n),rng.normal(size=(n,4))]),rng.normal(size=n))
    timings=[]
    for p in [0.,.3]:
        values={}
        for name,fn in [('original',original),('shortcut',fast)]:
            start=time.perf_counter()
            for _ in range(20):result=fn(cache,[.5,.7,p])
            values[name]=time.perf_counter()-start
        timings.append(dict(species_ratio=p,twenty_evaluations_seconds=values,speed_ratio=values['original']/values['shortcut']))
    record=dict(status='passed_exact_zero_species_gradient_shortcut_checks',cases=cases,derivatives=cases*3,
        maximum_dense_gradient_error=maximum,implementation_timings=timings,
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_reml_gradient_fast.py'),Path('scripts/matched_reml_gradient.py'),Path('scripts/cached_matched_likelihood.py')]},
        scope='Exact p=0 algebra; no near-zero cutoff, changed parameter bounds or omitted candidates. Fixed-point benchmark only; full optimization equivalence and timing remain pending. Running fits unchanged.')
    write_json(Path('metadata/matched_reml_zero_shortcut_checks_20260928.json'),record);print(record)


if __name__=='__main__':main()
