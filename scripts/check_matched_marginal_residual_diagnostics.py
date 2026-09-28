"""Dense projection oracle for residual variances, coefficients and moments."""
import json
from pathlib import Path
import numpy as np
from ancestral_chain_attempt import sha,write_json
from matched_marginal_residual_diagnostics import diagnostics
rng=np.random.default_rng(20260928);n=48
bg=np.repeat(np.arange(16),3);fam=bg//4;factor=rng.normal(size=(n,4))/2
x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n)
zb=np.eye(16)[bg];zf=np.eye(4)[fam]
results=[]
for ratios in [(0,0,0),(.3,0,.5),(0,.7,0),(.3,.7,.5)]:
    for scale in [.1,1.3,10.]:
        v=scale*np.r_[1.,ratios];V=v[0]*np.eye(n)+v[1]*zb@zb.T+v[2]*zf@zf.T+v[3]*factor@factor.T
        inv=np.linalg.inv(V);phi=np.linalg.inv(x.T@inv@x);hat=x@phi@x.T@inv
        covariance=(np.eye(n)-hat)@V@(np.eye(n)-hat).T
        beta=phi@x.T@inv@y;z=(y-x@beta)/np.sqrt(np.diag(covariance))
        r=diagnostics(bg,fam,factor,x,y,v)
        assert r['status']=='descriptive_marginal_residual_diagnostics'
        np.testing.assert_allclose(r['beta'],beta,rtol=1e-10,atol=1e-11)
        for name,power in [('mean',1),('second_raw_moment',2),('third_raw_moment',3),('fourth_raw_moment',4)]:
            np.testing.assert_allclose(r['summaries'][name],np.mean(z**power),rtol=1e-10,atol=1e-11)
        np.testing.assert_allclose(r['summaries']['quantiles'],np.quantile(z,r['summaries']['quantile_probabilities']),rtol=1e-10,atol=1e-11)
        assert r['summaries']['absolute_above_2']==sum(abs(z)>2)
        assert r['summaries']['absolute_above_3']==sum(abs(z)>3)
        results.append(dict(variances=v.tolist(),maximum_projection_diagonal_error=float(np.max(abs(np.diag(covariance)-np.diag(V-x@phi@x.T))))))
result=dict(status='passed_dense_residual_projection_checks',cases=results,
    pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_marginal_residual_diagnostics.py'),Path('scripts/matched_mixed_covariance.py')]},
    scope='Numerical residual projection checks, not model adequacy or calibrated tests.')
write_json(Path('metadata/matched_marginal_residual_checks_20260928.json'),result)
print(result['status'],len(results))
