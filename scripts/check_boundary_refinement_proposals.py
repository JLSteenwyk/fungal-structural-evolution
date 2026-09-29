#!/usr/bin/env python3
"""Direct-likelihood checks of saved boundary proposals, without gradient helper."""
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import profiled_ml
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    source=Path('metadata/boundary_refinement_diagnostic_completed_20260929.json')
    saved=json.loads(source.read_text())
    pp=Path('metadata/boundary_refinement_diagnostic_plan_20260929.json')
    plan=json.loads(pp.read_text())
    assert saved['plan_sha256']==sha(pp)
    bindings={**plan['pins'],str(source):sha(source),str(pp):sha(pp)}
    for path,digest in bindings.items():
        assert sha(path)==digest,path
    original=json.loads(Path(plan['fit']).read_text())['payload']
    with np.load(plan['inputs']) as a:
        numeric,bg,family,indices=a['matrix'],a['background'],a['family'],a['pattern_rows']
    with np.load(plan['factor']) as a:
        factor=a['factor'][indices]
    active=np.ptp(numeric[:,1:],axis=0)>1e-12
    scales=np.std(numeric[:,1:][:,active],axis=0)
    np.testing.assert_array_equal(active,original['active_covariates'])
    np.testing.assert_array_equal(scales,original['covariate_scales'])
    x=np.column_stack([np.ones(len(numeric)),numeric[:,1:][:,active]/scales]);y=numeric[:,0]
    evaluations=0
    def objective(theta):
        nonlocal evaluations
        evaluations+=1
        return profiled_ml(MatchedCovariance(bg,family,factor,1.,*np.expm1(theta)),x,y)['negative_profiled_ml']
    theta=np.asarray(original['log1p_ratios']);upper=np.log1p(original['maximum_ratio'])
    base=objective(theta)
    np.testing.assert_allclose(base,saved['original_objective'],rtol=0,atol=1e-8)
    assert len(saved['proposals'])==2
    results=[]
    for step,trial in zip([1e-5,1e-6],saved['proposals']):
        assert trial['hessian_step']==step
        point=np.asarray(trial['theta'])
        assert point[1]==0 and np.all(point[[0,2]]>1e-6) and np.all(point[[0,2]]<upper-1e-6)
        value=objective(point)
        np.testing.assert_allclose(value,trial['objective'],rtol=0,atol=1e-9)
        assert value<=base
        checks=[]
        for h in (1e-6,1e-7):
            g=np.empty(3)
            for axis in (0,2):
                plus=point.copy();minus=point.copy();plus[axis]+=h;minus[axis]-=h
                g[axis]=(objective(plus)-objective(minus))/(2*h)
            once=point.copy();twice=point.copy();once[1]=h;twice[1]=2*h
            f1,f2=objective(once),objective(twice)
            g[1]=((-3*value)+(4*f1)-f2)/(2*h)
            record=next(r for r in trial['finite_differences'] if r['step']==h)
            np.testing.assert_allclose(g,record['gradient'],rtol=1e-9,atol=1e-9)
            assert g[1]>0 and np.max(np.abs(g[[0,2]]))<=1e-3
            checks.append(dict(step=h,direct_gradient=g.tolist(),free_gradient_norm=float(np.max(np.abs(g[[0,2]]))),boundary_one_sided_gradient=float(g[1])))
        results.append(dict(hessian_step=step,theta=point.tolist(),objective=value,checks=checks))
    for path,digest in bindings.items():
        assert sha(path)==digest,path
    proof=dict(status='passed_saved_boundary_proposal_direct_checks',records=len(y),proposals=2,direct_evaluations=evaluations,checks=results,source_hashes=bindings,checker_sha256=sha(__file__),scope='Both saved proposal objectives and direct gradients checked using reconstructed design. No analytic-gradient or proposal helper called; direct likelihood library shared. Hessian construction not independently audited here. No fit replacement, global optimum or calibrated inference claim.')
    with Path('metadata/boundary_refinement_proposals_checked_20260929.json').open('x') as h:
        json.dump(proof,h,indent=2);h.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k not in ['checks','source_hashes']}),flush=True)


if __name__=='__main__':
    main()
