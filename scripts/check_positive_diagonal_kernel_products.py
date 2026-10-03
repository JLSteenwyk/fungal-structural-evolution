#!/usr/bin/env python3
"""Qualify general-diagonal raw/REML Grams against dense error contrasts."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse, linalg

from ancestral_chain_attempt import sha
from check_retained_shared_entity_candidates import inputs
from covariance_basis_context import ComponentKernelProducts
from independent_positive_diagonal_kernel_products import PositiveDiagonalKernelProducts
from reduced_covariance_basis import independent_diagnostics


def rejected(action):
    try: action()
    except (ValueError, AssertionError): return
    raise AssertionError('Malformed positive-diagonal audit inputs accepted')


def dense(incidence, factor, d, x):
    matrices = [np.diag(d), *[(z @ z.T).toarray() for z in incidence.values()], factor @ factor.T]
    # Complete QR gives an explicit orthonormal basis of error contrasts.
    q = linalg.qr(x, mode='full')[0]; contrasts = q[:,x.shape[1]:]
    residual = [contrasts.T @ k @ contrasts for k in matrices]
    raw = np.asarray([[np.sum(a*b) for b in matrices] for a in matrices])
    projected = np.asarray([[np.sum(a*b) for b in residual] for a in residual])
    return raw, projected


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(exist_ok=False);assert not args.receipt.exists()
    rng=np.random.default_rng(202610035);rows=[];cases=permutations=invalid=0
    maximum_raw=maximum_projected=0.;maximum_raw_bound_fraction=maximum_projected_bound_fraction=0.
    status=Counter();uniform_review=0;collapsed_kernel_review=0;transformed=0
    for exception,mode in itertools.product([False,True],['signed','unsigned']):
        source,selected,x,y,ops,audit,old,parent=inputs(exception,mode)
        labels=source['labels'];factor=source['factors'][audit['tree']];n=len(labels)
        target=sparse.eye(n,format='csr')
        operators={'target_node':target,'background_node':ops['background_node']}
        if exception:operators['model_pair']=ops['model_pair']
        operators['family_intercept']=ops['family_intercept']
        patterns={'uniform_one':np.ones(n),'constant_two':np.full(n,2.),
            'alternating':np.where(np.arange(n)%2,.5,2.),'positive_ramp':np.linspace(.125,3.,n),
            'near_uniform':1.+np.arange(n)*2.**-40,'wide_diagonal':np.geomspace(1e-4,1e4,n)}
        for pattern,d in patterns.items():
            primary=ComponentKernelProducts(labels,operators,d).tree(factor).audit(x)
            independent=PositiveDiagonalKernelProducts(labels,operators,d).tree(factor).project(x)
            reference=dense(operators,factor,d,x)
            raw=np.asarray(primary['raw_gram']);projected=np.asarray(primary['projected_gram'])
            raw_error=np.asarray(primary['raw_roundoff_envelope']);projected_error=np.asarray(primary['projected_roundoff_envelope'])
            for kind,a,b,e,f in [('raw',raw,reference[0],raw_error,independent['raw_error']),
                               ('projected',projected,reference[1],projected_error,independent['projected_error'])]:
                error=float(np.max(abs(a-b)));bound=e+f
                assert np.all(abs(a-b)<=bound+64*np.finfo(float).tiny),(kind,pattern,error)
                assert np.all(abs(a-independent[kind])<=bound+64*np.finfo(float).tiny),(kind,pattern)
                np.testing.assert_allclose(a,b,rtol=3e-9,atol=2e-8)
                np.testing.assert_allclose(a,independent[kind],rtol=3e-9,atol=2e-8)
                if kind=='raw':
                    maximum_raw=max(maximum_raw,error)
                    maximum_raw_bound_fraction=max(maximum_raw_bound_fraction,float(np.max(np.divide(abs(a-b),bound,out=np.zeros_like(bound),where=bound>0))))
                else:
                    maximum_projected=max(maximum_projected,error)
                    maximum_projected_bound_fraction=max(maximum_projected_bound_fraction,float(np.max(np.divide(abs(a-b),bound,out=np.zeros_like(bound),where=bound>0))))
            # The envelope and rank policy are reconstructed, not copied from uniform audits.
            for pfield,kind,ekind in [('raw_diagnostics','raw','raw_error'),('reml_diagnostics','projected','projected_error')]:
                check=independent_diagnostics(independent[kind],independent[ekind],primary['kernel_names'])
                assert primary[pfield]['disposition']==check['disposition'],(pattern,pfield)
                assert primary[pfield]['rank']==check['rank']
                status[primary[pfield]['disposition']]+=1
            if pattern in ['uniform_one','constant_two']:
                assert primary['raw_diagnostics']['disposition']=='dependent_covariance_bases_require_review'
                assert primary['reml_diagnostics']['disposition']=='dependent_covariance_bases_require_review'
                uniform_review+=1
            permutation=rng.permutation(n)
            perm=PositiveDiagonalKernelProducts(labels[permutation],{k:z[permutation] for k,z in operators.items()},d[permutation]).tree(factor[permutation]).project(x[permutation])
            for kind in ['raw','projected']:
                np.testing.assert_allclose(perm[kind],independent[kind],rtol=3e-9,atol=2e-8)
            permutations+=1
            # Full-rank column operations preserve the REML error-contrast subspace.
            transform=rng.normal(size=(x.shape[1],x.shape[1]));transform+=4*np.eye(x.shape[1])
            changed=PositiveDiagonalKernelProducts(labels,operators,d).tree(factor).project(x@transform)
            np.testing.assert_allclose(changed['projected'],independent['projected'],rtol=3e-9,atol=2e-8);transformed+=1
            rows.append(dict(pair_exception=exception,loading_mode=mode,pattern=pattern,
                primary=primary,independent={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in independent.items()},
                dense_raw=reference[0].tolist(),dense_projected=reference[1].tolist(),scientific_eligibility=False))
            cases+=1
        for d in [np.zeros(n),-np.ones(n),np.full(n,np.nan),np.full(n,np.inf),np.ones(n+1)]:
            rejected(lambda:PositiveDiagonalKernelProducts(labels,operators,d));invalid+=1
        bank=PositiveDiagonalKernelProducts(labels,operators,np.linspace(.5,2.,n))
        for design in [x[:,:1]*0,np.column_stack([x,x[:,0]]),np.ones((n,n)),x[:-1],np.full_like(x,np.nan)]:
            rejected(lambda:bank.tree(factor).project(design));invalid+=1
        for f in [np.full_like(factor,np.nan),factor[:-1],np.ones(n)]:
            rejected(lambda:bank.tree(f));invalid+=1
        bad=sparse.csr_matrix(np.ones((n,1)))
        rejected(lambda:PositiveDiagonalKernelProducts(labels,{'cross_component':bad},np.ones(n)));invalid+=1
        rejected(lambda:PositiveDiagonalKernelProducts(labels,{'residual':target},np.ones(n)));invalid+=1
        # A zero retained entity stays explicit and requires a norm review.
        zero={**operators,'explicit_zero':sparse.csr_matrix((n,0))}
        primary=ComponentKernelProducts(labels,zero,np.linspace(.5,2.,n)).tree(factor).audit(x)
        assert primary['raw_diagnostics']['disposition']=='unresolved_kernel_norm_requires_review'
        independent=PositiveDiagonalKernelProducts(labels,zero,np.linspace(.5,2.,n)).tree(factor).project(x)
        np.testing.assert_allclose(independent['raw'],primary['raw_gram'],rtol=3e-9,atol=2e-8)
        collapsed_kernel_review+=1
    assert (cases,permutations,transformed,uniform_review,collapsed_kernel_review,invalid)==(24,24,24,8,4,60)
    artifact=args.output/'positive_diagonal_gram_cases.json'
    artifact.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    modules=['check_positive_diagonal_kernel_products','independent_positive_diagonal_kernel_products',
        'covariance_basis_context','covariance_basis_audit','reduced_covariance_basis',
        'check_retained_shared_entity_candidates','covariance_exact_folds_v2','run_positive_diagonal_kernel_software_stage']
    result=dict(status='passed_positive_diagonal_latent_raw_reml_gram_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),synthetic_diagonal_cases=cases,
        independent_row_permutation_cases=permutations,full_rank_design_transform_cases=transformed,
        constant_diagonal_dependency_reviews=uniform_review,explicit_zero_kernel_reviews=collapsed_kernel_review,
        invalid_input_cases_rejected=invalid,diagnostic_status_counts=dict(status),
        maximum_dense_raw_gram_error=maximum_raw,maximum_dense_projected_gram_error=maximum_projected,
        maximum_raw_error_fraction_of_combined_roundoff_bound=maximum_raw_bound_fraction,
        maximum_projected_error_fraction_of_combined_roundoff_bound=maximum_projected_bound_fraction,
        unchanged_comparison_rtol=3e-9,unchanged_comparison_atol=2e-8,
        source_hashes={'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules},
        artifacts={str(artifact):sha(artifact)},scientific_eligibility=False,
        raw_reml_basis_qualification_complete=False,
        scope='Complete declared24positive diagonal q5/q6 signed/unsigned synthetic Gram checks against component arithmetic, independent latent overlaps and explicit dense QR error contrasts. Fresh bounds and rank reviews preserved; no uniform envelope inherited, basis deletion, real full-source weighted qualification, variance fit or biological acceptance.')
    with args.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','scope']}),flush=True)


if __name__=='__main__':main()
