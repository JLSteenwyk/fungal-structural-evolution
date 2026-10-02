#!/usr/bin/env python3
"""Compare cached block and independent latent-space residual kernel products."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from check_covariance_basis_audit import dense_check
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(553);errors=[];n=24;labels=np.repeat(np.arange(4),6)
    for trial in range(24):
        operators={}
        for name,width in [('gene',12),('model',8),('family',4)]:
            z=np.zeros((n,width));size=width//4
            for block in range(4):z[block*6:(block+1)*6,block*size:(block+1)*size]=rng.integers(-2,3,size=(6,size))/2
            # Sparse empty columns must not alter numerics or become dense.
            operators[name]=sparse.hstack([sparse.csr_matrix(z),sparse.csr_matrix((n,1000))],format='csr')
        factor=rng.normal(size=(n,trial%5));x=np.column_stack([np.ones(n),rng.normal(size=(n,2))])
        if trial%2:x=x*np.array([1e200,1e-200,1.])
        bank=ComponentKernelProducts(labels,operators,np.ones(n))
        independent=IndependentKernelProducts(labels,operators)
        for multiplier in [1.,1.25]:
            f=factor*multiplier;actual=bank.tree(f).audit(x);expected=independent.tree(f).project(x)
            oracle,_,_,_=dense_check(labels,operators,f,np.ones(n),x)
            for key,reference in [('raw_gram',expected['raw']),('projected_gram',expected['projected'])]:
                np.testing.assert_allclose(actual[key],reference,rtol=2e-10,atol=2e-9)
                np.testing.assert_allclose(actual[key],oracle[key],rtol=2e-10,atol=2e-9)
            np.testing.assert_allclose(actual['projected_roundoff_envelope'],expected['projected_error'],rtol=2e-9,atol=1e-20)
            errors.append(float(np.linalg.norm(np.asarray(actual['projected_gram'])-expected['projected'])/
                max(np.linalg.norm(expected['projected']),np.finfo(float).tiny)))
            # Changing designs must leave reusable raw kernel products intact.
            before=bank.tree(f).raw.copy();bank.tree(f).audit(np.ones((n,1)))
            np.testing.assert_array_equal(before,bank.tree(f).raw)
    invalid=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,ArithmeticError):invalid.append(name)
        else:raise AssertionError('Invalid context accepted: '+name)
    reject('cross_component_entity',lambda:ComponentKernelProducts(labels,{'z':sparse.csr_matrix(np.ones((n,1)))},np.ones(n)))
    reject('nonpositive_residual',lambda:ComponentKernelProducts(labels,operators,np.zeros(n)))
    reject('nonfinite_entity',lambda:ComponentKernelProducts(labels,{'z':sparse.csr_matrix(np.full((n,1),np.nan))},np.ones(n)))
    reject('reserved_entity_name',lambda:ComponentKernelProducts(labels,{'residual':sparse.eye(n)},np.ones(n)))
    reject('nonfinite_factor',lambda:bank.tree(np.full((n,1),np.inf)))
    reject('rank_deficient_design',lambda:bank.tree(factor).audit(np.ones((n,2))))
    reject('no_residual_dimensions',lambda:bank.tree(factor).audit(np.eye(n)))
    sources=[Path(__file__),*[Path('scripts')/v for v in ['covariance_basis_context.py','covariance_basis_independent.py',
        'covariance_basis_audit.py','check_covariance_basis_audit.py']]]
    result=dict(status='passed_reusable_covariance_kernel_context_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(),independent_dense_and_latent_cases=len(errors),
        maximum_projected_gram_relative_frobenius_error=max(errors),
        empty_global_columns_removed_without_kernel_change=True,raw_cache_unchanged_by_design=True,
        rejected_invalid_inputs=invalid,source_hashes={str(p):sha(p) for p in sources},
        scope='Synthetic numerical contracts only. Block-space, latent-space and dense-null-space '
              'algebra agree under signed loadings, rank-zero species factors, reuse and extreme '
              'design units. Not production-cohort qualification, a biological pilot or a variance fit.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
