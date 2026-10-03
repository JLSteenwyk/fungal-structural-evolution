#!/usr/bin/env python3
"""Dense/high-precision/locality and explicit precision-fallback contracts."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
from scipy import sparse
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood,audit_local_curvature
from independent_shared_entity_likelihood_fast import ComponentSpectralLikelihood
from independent_shared_entity_optimizer import compare_independent_optimizer
from shared_entity_likelihood import SharedEntityLikelihood
from fit_shared_entity_likelihood import fit_shared_entity
from check_shared_entity_likelihood import dense
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    rng=np.random.default_rng(753);cases=0;objective_error=gradient_error=0.
    keys=['negative_profiled_likelihood','beta','conditional_beta_covariance','profiled_scale',
        'gradient','residual_quadratic','kernel_traces_in_score','residual_kernel_energies']
    for trial in range(24):
        n=18+trial%13;labels=np.arange(n)%3;z={}
        for name in ['gene','model','family']:
            values=np.zeros((n,12))
            for column in range(9):values[labels==column%3,column]=rng.normal(size=(labels==column%3).sum())
            z[name]=sparse.csr_matrix(values)
        if trial%4==0:z['cancelled']=sparse.csr_matrix((n,17))
        f=rng.normal(size=(n,trial%5));d=rng.uniform(.4,1.5,n);x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n)
        ratios=rng.uniform(.01,3.,len(z)+1)
        if trial%2==0:ratios[0]=0.
        if trial%6==0:ratios[:]=0.
        slow=IndependentSharedEntityLikelihood(labels,z,f,d,x,y,column_batch=2)
        fast=ComponentSpectralLikelihood(labels,z,f,d,x,y,column_batch=2)
        for method in ['ml','reml']:
            expected=slow.evaluate(ratios,method);actual=fast.evaluate(ratios,method);reference=dense(labels,z,f,d,x,y,ratios,method)
            for key in keys:np.testing.assert_allclose(actual[key],expected[key],rtol=2e-8,atol=2e-9)
            for key in reference:np.testing.assert_allclose(actual[key],reference[key],rtol=2e-8,atol=2e-9)
            assert actual['largest_dense_kernel_column_batch']<=2 and actual['scientific_eligibility'] is False
            objective_error=max(objective_error,abs(actual['negative_profiled_likelihood']-reference['negative_profiled_likelihood']))
            gradient_error=max(gradient_error,float(np.max(abs(actual['gradient']-reference['gradient']))))
            if trial<4:
                for j,value in enumerate(ratios):
                    h=1e-4*max(1.,value)
                    def point(distance):
                        changed=ratios.copy();changed[j]+=distance*h
                        return fast.evaluate(changed,method)['negative_profiled_likelihood']
                    numeric=(-25*point(0)+48*point(1)-36*point(2)+16*point(3)-3*point(4))/(12*h) if value<2*h else (point(-2)-8*point(-1)+8*point(1)-point(2))/(12*h)
                    np.testing.assert_allclose(actual['gradient'][j],numeric,rtol=3e-5,atol=3e-7)
            permutation=rng.permutation(n)
            permuted=ComponentSpectralLikelihood(labels[permutation],{k:v[permutation] for k,v in z.items()},f[permutation],d[permutation],x[permutation],y[permutation],column_batch=3).evaluate(ratios,method)
            for key in keys:np.testing.assert_allclose(actual[key],permuted[key],rtol=2e-8,atol=2e-9)
            cases+=1
    # Exact annihilation of a kernel by the intercept must use explicit norms.
    n=18;x=np.column_stack([np.ones(n),rng.normal(size=n)]);y=rng.normal(size=n)
    annihilated=ComponentSpectralLikelihood(np.zeros(n),{'intercept':sparse.csr_matrix(np.ones((n,1)))},np.empty((n,0)),np.ones(n),x,y).evaluate([10.,0.])
    assert annihilated['explicit_global_whitening_fallback_batches']==1
    assert annihilated['kernel_traces_in_score'][0]>=0
    np.testing.assert_allclose(annihilated['gradient'][0],0.,atol=1e-10)
    # Strong signed high-precision reference, independent of both float paths.
    z=np.array([[1,0],[1,0],[-1,0],[0,1],[0,1],[0,-1]],dtype=float)
    f=np.array([[.1,.2],[.2,-.1],[.3,.3],[.4,-.2],[.5,.5],[.6,.4]])
    d=np.array([.2,.5,.7,.3,.6,.4]);x=np.column_stack([np.ones(6),np.arange(1,7)/10]);y=np.array([.3,-.8,.1,.6,-.3,.2])
    def matrix(values):return mp.matrix([[mp.mpf(str(v)) for v in row] for row in values])
    precision={}
    with mp.workdps(80):
        mz=matrix(z);mf=matrix(f);mx=matrix(x);my=mp.matrix([mp.mpf(str(v)) for v in y])
        kernels=[mz*mz.T,mf*mf.T];v=mp.diag([mp.mpf(str(v)) for v in d])+10000*kernels[0]+3000*kernels[1]
        iv=v**-1;info=mx.T*iv*mx;ii=info**-1;beta=ii*mx.T*iv*my;residual=my-mx*beta;alpha=iv*residual
        q=(residual.T*alpha)[0];projection=iv-iv*mx*ii*mx.T*iv
        for method in ['ml','reml']:
            df=4 if method=='reml' else 6
            objective=(mp.log(mp.det(v))+(mp.log(mp.det(info)) if method=='reml' else 0)+df*(1+mp.log(2*mp.pi*q/df)))/2
            gradient=np.array([float((sum(((projection if method=='reml' else iv)*k)[i,i] for i in range(6))-df*(alpha.T*k*alpha)[0]/q)/2) for k in kernels])
            actual=ComponentSpectralLikelihood(['A']*3+['B']*3,{'signed':sparse.csr_matrix(z)},f,d,x,y,column_batch=1).evaluate([10000.,3000.],method)
            np.testing.assert_allclose(actual['negative_profiled_likelihood'],float(objective),rtol=2e-9,atol=2e-9)
            np.testing.assert_allclose(actual['gradient'],gradient,rtol=2e-6,atol=2e-10)
            precision[method]=dict(objective_absolute_error=abs(actual['negative_profiled_likelihood']-float(objective)),
                maximum_gradient_absolute_error=float(np.max(abs(actual['gradient']-gradient))))
    # Resource-bounded software benchmark: 6,000 synthetic rows, 300 blocks,
    # three sparse incidence arrays and rank-eight factor, no biological data.
    n=6000;block_size=20;block_count=n//block_size;labels=np.repeat(np.arange(block_count),block_size);ops={}
    row_indices=np.repeat(np.arange(n),12)
    column_indices=np.repeat(np.repeat(np.arange(block_count)*12,block_size),12)+np.tile(np.arange(12),n)
    for name in ['gene','model','pair']:ops[name]=sparse.csr_matrix((rng.normal(size=len(row_indices)),(row_indices,column_indices)),shape=(n,block_count*12))
    f=rng.normal(size=(n,8));d=np.ones(n);x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n);ratios=np.array([.3,.7,1.,.2])
    fast=ComponentSpectralLikelihood(labels,ops,f,d,x,y);slow=IndependentSharedEntityLikelihood(labels,ops,f,d,x,y)
    started=time.perf_counter();a=fast.evaluate(ratios);fast_seconds=time.perf_counter()-started
    started=time.perf_counter();b=slow.evaluate(ratios);slow_seconds=time.perf_counter()-started
    for key in keys:np.testing.assert_allclose(a[key],b[key],rtol=2e-8,atol=2e-8)
    conceptual_global_cells=n*sum(z.shape[1] for z in slow.incidence.values())
    assert a['component_local_entity_cells']*100<conceptual_global_cells and a['explicit_global_whitening_fallback_batches']==0
    # The existing numerical candidate/curvature protocol accepts the subtype.
    rng=np.random.default_rng(473);n=36;labels=np.repeat(np.arange(3),12);ops={}
    for name in ['gene','family']:
        z=np.zeros((n,9))
        for column in range(9):z[labels==column%3,column]=rng.normal(size=12)
        ops[name]=sparse.csr_matrix(z)
    f=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)]);y=ops['gene']@rng.normal(size=9)+.3*(ops['family']@rng.normal(size=9))+rng.normal(size=n)
    fitted=fit_shared_entity(SharedEntityLikelihood(labels,ops,f,np.ones(n),x,y),'reml',gradient_tolerance=3e-6)
    fast=ComponentSpectralLikelihood(labels,ops,f,np.ones(n),x,y,column_batch=3)
    assert audit_local_curvature(fast,fitted,gradient_atol=3e-6)['status']=='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits'
    assert compare_independent_optimizer(fast,fitted,gradient_tolerance=3e-6)['status']=='independent_multistart_searches_agree_pending_inferential_calibration'
    receipt=dict(status='passed_component_local_spectral_likelihood_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        dense_and_original_streamed_ml_reml_cases=cases,maximum_dense_objective_absolute_error=objective_error,
        maximum_dense_gradient_absolute_error=gradient_error,precision_digits=80,strong_signed_precision_errors=precision,
        explicit_precision_fallback_exercised=True,row_permutations_and_bounded_batches_passed=True,
        independent_curvature_and_optimizer_subtype_path_passed=True,
        synthetic_locality_benchmark=dict(rows=6000,components=300,factor_rank=8,conceptual_global_entity_whitening_cells=conceptual_global_cells,
            actual_component_local_entity_cells=a['component_local_entity_cells'],global_fallback_batches=a['explicit_global_whitening_fallback_batches'],
            cached_component_kernel_bytes=a['cached_component_kernel_bytes'],component_evaluation_seconds=fast_seconds,original_streamed_evaluation_seconds=slow_seconds,
            speed_ratio_at_observation=slow_seconds/fast_seconds),
        source_hashes={p:sha(p) for p in ['scripts/check_component_spectral_likelihood.py','scripts/independent_shared_entity_likelihood_fast.py',
            'scripts/independent_shared_entity_likelihood.py','scripts/independent_shared_entity_optimizer.py','scripts/shared_entity_likelihood.py',
            'scripts/fit_shared_entity_likelihood.py','scripts/check_shared_entity_likelihood.py','scripts/shared_entity_covariance.py']},
        scope='Synthetic numerical contracts and a resource-bounded software locality benchmark only. '
            'No biological pilot, production cohort timing, full-stage runtime estimate or calibrated inference. '
            'Component-local trace differences use a conservative cancellation screen and original explicit '
            'global whitening/norm fallback; no jitter or eigenvalue/trace clipping.')
    with args.output.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
