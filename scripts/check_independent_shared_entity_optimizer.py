#!/usr/bin/env python3
"""Independent SLSQP contracts with dense objectives and explicit review states."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from shared_entity_likelihood import SharedEntityLikelihood
from fit_shared_entity_likelihood import fit_shared_entity
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood
from independent_shared_entity_optimizer import compare_independent_optimizer
from check_shared_entity_likelihood import dense
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(473);cases=passed=0;maximum_error=0.;statuses=[]
    for trial in range(3):
        n=36;labels=np.repeat(np.arange(3),12);z={}
        for name in ['gene','family']:
            values=np.zeros((n,9))
            for column in range(9):values[labels==column%3,column]=rng.normal(size=12)
            z[name]=sparse.csr_matrix(values)
        f=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)]);d=np.ones(n)
        y=z['gene']@rng.normal(size=9)+.3*(z['family']@rng.normal(size=9))+rng.normal(size=n)
        original=SharedEntityLikelihood(labels,z,f,d,x,y);independent=IndependentSharedEntityLikelihood(labels,z,f,d,x,y,column_batch=3)
        for method in ['ml','reml']:
            candidate=fit_shared_entity(original,method,gradient_tolerance=3e-6)
            result=compare_independent_optimizer(independent,candidate,gradient_tolerance=3e-6)
            assert len(result['starts'])==3,result
            for start in result['starts']:
                if 'final_coordinates' in start:
                    ratios=np.expm1(start['final_coordinates'])/np.asarray(candidate['kernel_normalization'])
                    reference=dense(labels,z,f,d,x,y,ratios,method)['negative_profiled_likelihood']
                    error=abs(reference-start['objective']);maximum_error=max(maximum_error,error)
                    assert error<1e-7
            passed+=result['status']=='independent_multistart_searches_agree_pending_inferential_calibration'
            statuses.append(result['status']);cases+=1
            if trial==0:
                assert result['status']=='independent_multistart_searches_agree_pending_inferential_calibration',result
                exhausted=compare_independent_optimizer(independent,candidate,max_iterations=1,max_evaluations=1,gradient_tolerance=3e-6)
                assert exhausted['status']=='independent_shared_entity_search_requires_review'
                assert all(s['disposition']!='independent_search_agrees_pending_calibration' and s['search_evaluations']<=1 for s in exhausted['starts'])
                altered=copy.deepcopy(candidate);altered['negative_profiled_likelihood']+=1.
                rejected=compare_independent_optimizer(independent,altered)
                assert rejected['status']=='independent_shared_entity_search_requires_review' and not rejected['starts']
    assert passed>=2
    sources=['scripts/check_independent_shared_entity_optimizer.py','scripts/independent_shared_entity_optimizer.py',
        'scripts/independent_shared_entity_likelihood.py','scripts/shared_entity_likelihood.py',
        'scripts/shared_entity_covariance.py','scripts/fit_shared_entity_likelihood.py',
        'scripts/check_shared_entity_likelihood.py','scripts/covariance_basis_context.py','scripts/covariance_basis_audit.py']
    receipt=dict(status='passed_independent_spectral_multistart_optimizer_contracts',independent_ml_reml_cases=cases,
        candidates_with_all_independent_starts_agreeing=passed,candidate_statuses=statuses,maximum_dense_objective_error=maximum_error,
        independent_budget_exhaustion_retained_as_review=True,altered_candidate_refused_before_search=True,
        source_hashes={p:sha(p) for p in sources},scope='Synthetic complete numerical optimizer fixtures only. '
            'Dense objective checks of independent spectral SLSQP starts. Agreement is not a global '
            'optimality proof, production biological fit, pilot or calibrated uncertainty.')
    with a.output.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
