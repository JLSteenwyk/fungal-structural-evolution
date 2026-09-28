"""Compare batched REML to scalar and independently assembled dense covariance."""
import itertools
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from batched_matched_reml import evaluate
from matched_mixed_covariance import MatchedCovariance, profiled_reml


def main():
    threadpool_limits(1)
    rng = np.random.default_rng(20260928)
    n = 24; bg = np.repeat(np.arange(8),3); fam = bg//2
    x = np.column_stack([np.ones(n),rng.normal(size=n),rng.normal(size=n)])
    factor = rng.normal(size=(n,4))/2
    responses = rng.normal(size=(n,7))
    # Retain a failed zero response beside valid responses; no batch deletion.
    responses = np.column_stack([responses,np.zeros(n)])
    errors = []; comparisons = 0
    for rank in [0,4]:
        f = factor[:,:rank]
        for ratios in itertools.product([0.,2.],repeat=3):
            cov = MatchedCovariance(bg,fam,f,1.,*ratios)
            result = evaluate(cov,x,responses)
            assert result['valid'].tolist() == [True]*7+[False]
            assert np.isnan(result['negative_profiled_reml'][-1])
            v = np.eye(n)+ratios[0]*(bg[:,None]==bg[None,:])+ratios[1]*(fam[:,None]==fam[None,:])+ratios[2]*(f@f.T)
            invx = np.linalg.solve(v,x); info=x.T@invx
            for j in range(7):
                y = responses[:,j]
                scalar = profiled_reml(cov,x,y)
                beta = np.linalg.solve(info,invx.T@y)
                residual = y-x@beta
                scale = (residual@np.linalg.solve(v,residual))/(n-x.shape[1])
                dense = .5*(np.linalg.slogdet(v)[1]+np.linalg.slogdet(info)[1]+(n-x.shape[1])*(1+np.log(2*np.pi*scale)))
                np.testing.assert_allclose(result['negative_profiled_reml'][j],dense,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(result['negative_profiled_reml'][j],scalar['negative_profiled_reml'],rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(result['beta'][:,j],beta,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(result['conditional_beta_covariance'][j],scale*np.linalg.inv(info),rtol=1e-12,atol=1e-12)
                errors.append(abs(result['negative_profiled_reml'][j]-dense));comparisons+=1
    invalid=[(x,responses[:,0]),(x[:-1],responses),(np.ones((n,2)),responses),(x,np.full((n,1),np.nan))]
    for xx,yy in invalid:
        try:evaluate(cov,xx,yy)
        except ValueError:pass
        else:raise AssertionError('Invalid input accepted')
    # Small implementation benchmark only; it does not time independent optimizers.
    n=600;bg=np.repeat(np.arange(100),6);fam=bg//5
    x=np.column_stack([np.ones(n),rng.normal(size=(n,4))]);f=rng.normal(size=(n,20))/5
    y=rng.normal(size=(n,64));cov=MatchedCovariance(bg,fam,f,1.,.5,.7,.3)
    timings=[]
    for _ in range(3):
        start=time.perf_counter();batch=evaluate(cov,x,y);batched=time.perf_counter()-start
        start=time.perf_counter();scalar=[profiled_reml(cov,x,y[:,j]) for j in range(y.shape[1])];single=time.perf_counter()-start
        np.testing.assert_allclose(batch['negative_profiled_reml'],[r['negative_profiled_reml'] for r in scalar],rtol=1e-12,atol=1e-10)
        timings.append(dict(batched_seconds=batched,scalar_seconds=single,speed_ratio=single/batched))
    record=dict(status='passed_batched_scalar_and_dense_reml_checks',dense_response_comparisons=comparisons,
        invalid_input_cases=len(invalid),zero_response_retained_as_invalid=True,
        maximum_dense_objective_error=max(errors),benchmark=dict(observations=600,responses=64,design_columns=5,species_rank=20,trials=timings),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/batched_matched_reml.py'),Path('scripts/matched_mixed_covariance.py')]},
        scope='Fixed common covariance only. Timing excludes covariance construction and per-response ratio optimization; cannot extrapolate to full refit speed or calibrated uncertainty.')
    write_json(Path('metadata/batched_matched_reml_checks_20260928.json'),record)
    print(record)


if __name__=='__main__':main()
