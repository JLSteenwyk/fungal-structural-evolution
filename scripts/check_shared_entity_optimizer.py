#!/usr/bin/env python3
"""Independent dense optimizer and bounded-search disposition contracts."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import minimize
from shared_entity_likelihood import SharedEntityLikelihood
from fit_shared_entity_likelihood import fit_shared_entity
from check_shared_entity_likelihood import dense
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(473);checked=qualified=0;maximum_error=0.;statuses=[]
    for trial in range(3):
        n=36;labels=np.repeat(np.arange(3),12);operators={}
        for name in ['gene','family']:
            z=np.zeros((n,9))
            for column in range(9):z[labels==column%3,column]=rng.normal(size=12)
            operators[name]=sparse.csr_matrix(z)
        factor=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)]);d=np.ones(n)
        y=operators['gene']@rng.normal(size=9)+.3*(operators['family']@rng.normal(size=9))+rng.normal(size=n)
        likelihood=SharedEntityLikelihood(labels,operators,factor,d,x,y)
        for method in ['ml','reml']:
            result=fit_shared_entity(likelihood,method,gradient_tolerance=3e-6)
            assert 'variance_ratios' in result
            expected=dense(labels,operators,factor,d,x,y,np.asarray(result['variance_ratios']),method)
            np.testing.assert_allclose(result['negative_profiled_likelihood'],expected['negative_profiled_likelihood'],rtol=2e-10,atol=2e-9)
            norms=np.asarray(result['kernel_normalization']);upper=np.log1p(result['maximum_scaled_variance'])
            def independent(point):return dense(labels,operators,factor,d,x,y,np.expm1(point)/norms,method)['negative_profiled_likelihood']
            reference=min([minimize(independent,np.full(len(norms),v),method='Powell',bounds=[(0.,upper)]*len(norms),
                options=dict(xtol=1e-10,ftol=1e-12,maxiter=500)) for v in [0.,np.log(2.)]],key=lambda r:r.fun)
            error=abs(reference.fun-result['negative_profiled_likelihood'])
            assert error<1e-5,(trial,method,error,result['starts'])
            maximum_error=max(maximum_error,error);checked+=1
            statuses.append(result['status']);qualified+=result['status']=='optimized_shared_entity_candidate_pending_independent_audit'
            assert result['scientific_eligibility'] is False and result['coefficient_covariance_is_conditional']
            assert len(result['starts'])==3
            capped=fit_shared_entity(likelihood,method,maximum_scaled_variance=.001)
            assert capped['status']=='shared_entity_optimizer_candidate_requires_review'
            assert any(s.get('upper_boundary_indices') for s in capped['starts'])
    assert qualified>0,'No candidate exercised the fully converged multi-start path'
    # L-BFGS-B's displacement-based projection can stop immediately in a
    # very narrow box despite a large true boundary KKT violation.
    narrow=fit_shared_entity(likelihood,maximum_scaled_variance=1e-7)
    assert narrow['status']=='shared_entity_optimizer_candidate_requires_review'
    assert any(s.get('optimizer_success') and s.get('maximum_projected_gradient',0)>1e-6 for s in narrow['starts'])
    exhausted=fit_shared_entity(likelihood,max_iterations=1,max_evaluations=1)
    assert exhausted['status']=='all_shared_entity_optimizer_searches_require_review'
    assert all(s['disposition']=='failed_search_requires_review' for s in exhausted['starts'])
    identity=SharedEntityLikelihood(labels,{'target':sparse.eye(n,format='csr')},factor,d,x,y)
    try:fit_shared_entity(identity)
    except ValueError:pass
    else:raise AssertionError('Unidentified target/residual variance accepted')
    result=dict(status='passed_shared_entity_multistart_optimizer_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        independent_dense_powell_cases=checked,multistart_candidates_pending_audit=qualified,candidate_statuses=statuses,
        maximum_independent_optimizer_objective_error=maximum_error,upper_variance_caps_retained_as_review=True,
        exhausted_evaluation_budgets_retained_as_failure=True,target_residual_alias_rejected=True,
        narrow_box_false_optimizer_success_retained_as_review=True,
        source_hashes={str(p):sha(p) for p in [Path(__file__),Path('scripts/fit_shared_entity_likelihood.py'),
            Path('scripts/shared_entity_likelihood.py'),Path('scripts/check_shared_entity_likelihood.py'),
            Path('scripts/covariance_basis_context.py'),Path('scripts/covariance_basis_audit.py')]},
        scope='Synthetic software optimizer contracts only. Independent dense Powell solutions and '
              'explicit boundary/failure/identifiability tests. No production fit, biological pilot, '
              'independent full-data optimization acceptance or calibrated uncertainty.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
