#!/usr/bin/env python3
"""Fixed-face local proposals with explicit one-sided boundary derivatives."""
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/boundary_refinement_diagnostic_plan_20260929.json')
    plan = json.loads(pp.read_text())
    for path,digest in plan['pins'].items():
        assert sha(path) == digest,path
    saved = json.loads(Path(plan['fit']).read_text())['payload']
    with np.load(plan['inputs']) as a:
        values,bg,family,rows = a['matrix'],a['background'],a['family'],a['pattern_rows']
    with np.load(plan['factor']) as a:
        factor = a['factor'][rows]
    covariates = values[:,1:][:,saved['active_covariates']]
    x = np.column_stack((np.ones(len(values)),covariates/np.std(covariates,axis=0)))
    y=values[:,0]
    cache=CachedMatchedML(bg,family,factor,x,y)
    def objective(theta):
        return profiled_ml(MatchedCovariance(bg,family,factor,1.,*np.expm1(theta)),x,y)['negative_profiled_ml']
    def gradient(theta):
        return evaluate_gradient(cache,np.expm1(theta))['log1p_ratio_gradient']
    theta=np.asarray(saved['log1p_ratios']);upper=np.log1p(saved['maximum_ratio'])
    fixed=np.flatnonzero(theta==0);active=np.flatnonzero(theta>0)
    assert fixed.tolist()==[1] and active.tolist()==[0,2]
    assert np.all(theta[active]>1e-5) and np.all(theta[active]<upper-1e-5)
    base=objective(theta);g=gradient(theta)
    np.testing.assert_allclose(base,saved['negative_profiled_ml'],rtol=1e-12,atol=1e-8)
    np.testing.assert_allclose(g,saved['analytic_gradient'],rtol=1e-9,atol=1e-9)
    def direct_checks(point):
        f0=objective(point);out=[]
        for h in (1e-6,1e-7):
            estimates=[]
            for j in range(3):
                offset=np.eye(3)[j]*h
                if j in fixed:
                    derivative=(-3*f0+4*objective(point+offset)-objective(point+2*offset))/(2*h)
                else:
                    derivative=(objective(point+offset)-objective(point-offset))/(2*h)
                estimates.append(float(derivative))
            projected=np.asarray(estimates);projected[fixed]=np.minimum(projected[fixed],0)
            out.append(dict(step=h,gradient=estimates,projected_norm=float(np.max(abs(projected))),passes=bool(np.max(abs(projected))<=1e-3)))
        return out
    trials=[]
    for h in (1e-5,1e-6):
        raw=np.column_stack([(gradient(theta+np.eye(3)[j]*h)[active]-gradient(theta-np.eye(3)[j]*h)[active])/(2*h) for j in active])
        hessian=.5*(raw+raw.T);eigenvalues=np.linalg.eigvalsh(hessian)
        assert np.all(eigenvalues>0)
        point=theta.copy();point[active]-=np.linalg.solve(hessian,g[active])
        assert np.all(point[active]>1e-6) and np.all(point[active]<upper-1e-6) and np.all(point[fixed]==0)
        value=objective(point);derivative=gradient(point)
        projected=derivative.copy();projected[fixed]=np.minimum(projected[fixed],0)
        fd=direct_checks(point)
        checks=dict(objective_not_worsened=bool(value<=base),analytic_projected_gradient_pass=bool(np.max(abs(projected))<=1e-3),both_direct_projected_gradients_pass=all(d['passes'] for d in fd))
        trials.append(dict(hessian_step=h,hessian=hessian.tolist(),eigenvalues=eigenvalues.tolist(),theta=point.tolist(),objective=value,gradient=derivative.tolist(),projected_gradient=projected.tolist(),finite_differences=fd,checks=checks))
    for path,digest in plan['pins'].items():
        assert sha(path)==digest,path
    result=dict(status='complete_boundary_case_local_diagnostic',records=len(y),original_theta=theta.tolist(),original_objective=base,original_gradient=g.tolist(),original_direct_checks=direct_checks(theta),proposals=trials,source_hashes=plan['pins'],plan_sha256=sha(pp),script_sha256=sha(__file__),scope='Single fixed-face diagnostic; family variance remains exactly zero. Second-order forward differences check its one-sided gradient. No production replacement; separate readback and upstream full audit required. No global optimum or inferential claim.')
    with Path(plan['proof']).open('x') as handle:
        json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(status=result['status'],original_gradient=g.tolist(),proposal_checks=[t['checks'] for t in trials])),flush=True)


if __name__=='__main__':
    main()
