#!/usr/bin/env python3
"""Dense, production and altered-candidate contracts for spectral replay."""
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy import sparse
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood,audit_candidate,audit_local_curvature
from shared_entity_likelihood import SharedEntityLikelihood
from fit_shared_entity_likelihood import fit_shared_entity
from check_shared_entity_likelihood import dense
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(1952);cases=0;maximum_objective=maximum_gradient=maximum_inverse=0.
    for trial in range(24):
        n=15+trial%11;labels=np.arange(n)%3;operators={}
        for name in ['gene','model','family']:
            z=np.zeros((n,12))
            for column in range(9):z[labels==column%3,column]=rng.normal(size=(labels==column%3).sum())
            operators[name]=sparse.csr_matrix(z)
        if trial%4==0:operators['cancelled']=sparse.csr_matrix((n,30))
        factor=rng.normal(size=(n,trial%5));diagonal=rng.uniform(.3,2.,n)
        x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n)
        ratios=rng.uniform(.1,2.,len(operators)+1)
        if trial%2==0:ratios[0]=0.
        if trial%6==0:ratios[:]=0.
        for method in ['ml','reml']:
            independent=IndependentSharedEntityLikelihood(labels,operators,factor,diagonal,x,y,column_batch=2)
            actual=independent.evaluate(ratios,method)
            expected=dense(labels,operators,factor,diagonal,x,y,ratios,method)
            producer=SharedEntityLikelihood(labels,operators,factor,diagonal,x,y).evaluate(ratios,method)
            for key in expected:
                np.testing.assert_allclose(actual[key],expected[key],rtol=2e-8,atol=2e-9)
                np.testing.assert_allclose(actual[key],producer[key],rtol=2e-8,atol=2e-9)
            maximum_objective=max(maximum_objective,abs(actual['negative_profiled_likelihood']-expected['negative_profiled_likelihood']))
            maximum_gradient=max(maximum_gradient,float(np.max(abs(actual['gradient']-expected['gradient']))))
            maximum_inverse=max(maximum_inverse,actual['inverse_relative_residual'])
            assert actual['largest_dense_kernel_column_batch']<=2 and actual['scientific_eligibility'] is False
            alternate=IndependentSharedEntityLikelihood(labels,operators,factor,diagonal,x,y,column_batch=32).evaluate(ratios,method)
            np.testing.assert_allclose(actual['gradient'],alternate['gradient'],rtol=2e-9,atol=2e-10)
            permutation=rng.permutation(n)
            permuted=IndependentSharedEntityLikelihood(labels[permutation],{k:z[permutation] for k,z in operators.items()},
                factor[permutation],diagonal[permutation],x[permutation],y[permutation]).evaluate(ratios,method)
            for key in expected:np.testing.assert_allclose(actual[key],permuted[key],rtol=2e-8,atol=2e-9)
            cases+=1
    units=np.array([1e100,1e-100,1.])
    baseline=IndependentSharedEntityLikelihood(labels,operators,factor,diagonal,x,y).evaluate(ratios)
    scaled=IndependentSharedEntityLikelihood(labels,operators,factor,diagonal,x*units,y).evaluate(ratios)
    np.testing.assert_allclose(baseline['gradient'],scaled['gradient'],rtol=2e-8,atol=2e-9)
    np.testing.assert_allclose(baseline['beta'],scaled['beta']*units,rtol=2e-8,atol=2e-9)
    np.testing.assert_allclose(baseline['conditional_beta_covariance'],scaled['conditional_beta_covariance']*units[:,None]*units[None,:],rtol=2e-8,atol=2e-9)
    # Strong signed case also has an independent 80-digit check in a companion checker.
    z=np.array([[1,0],[1,0],[-1,0],[0,1],[0,1],[0,-1]],dtype=float)
    f=np.array([[.1,.2],[.2,-.1],[.3,.3],[.4,-.2],[.5,.5],[.6,.4]])
    d=np.array([.2,.5,.7,.3,.6,.4]);sx=np.column_stack([np.ones(6),np.arange(1,7)/10]);sy=np.array([.3,-.8,.1,.6,-.3,.2])
    precision_errors={}
    def mp_matrix(values):return mp.matrix([[mp.mpf(str(v)) for v in row] for row in values])
    with mp.workdps(80):
        mz=mp_matrix(z);mf=mp_matrix(f);mx=mp_matrix(sx);my=mp.matrix([mp.mpf(str(v)) for v in sy])
        kernels=[mz*mz.T,mf*mf.T];v=mp.diag([mp.mpf(str(v)) for v in d])+10000*kernels[0]+3000*kernels[1]
        iv=v**-1;info=mx.T*iv*mx;ii=info**-1;mbeta=ii*mx.T*iv*my
        residual=my-mx*mbeta;alpha=iv*residual;quadratic=(residual.T*alpha)[0]
        projection=iv-iv*mx*ii*mx.T*iv
        for method in ['ml','reml']:
            df=4 if method=='reml' else 6
            objective=(mp.log(mp.det(v))+(mp.log(mp.det(info)) if method=='reml' else 0)+df*(1+mp.log(2*mp.pi*quadratic/df)))/2
            gradient=np.array([float((sum(((projection if method=='reml' else iv)*k)[i,i] for i in range(6))-df*(alpha.T*k*alpha)[0]/quadratic)/2) for k in kernels])
            actual=IndependentSharedEntityLikelihood(['A']*3+['B']*3,{'signed':sparse.csr_matrix(z)},f,d,sx,sy).evaluate([10000.,3000.],method)
            np.testing.assert_allclose(actual['negative_profiled_likelihood'],float(objective),rtol=2e-9,atol=2e-9)
            np.testing.assert_allclose(actual['gradient'],gradient,rtol=2e-6,atol=2e-10)
            precision_errors[method]=dict(objective_absolute_error=abs(actual['negative_profiled_likelihood']-float(objective)),
                maximum_gradient_absolute_error=float(np.max(abs(actual['gradient']-gradient))))
    for method in ['ml','reml']:
        spectral=IndependentSharedEntityLikelihood(['A']*3+['B']*3,{'signed':sparse.csr_matrix(z)},f,d,sx,sy).evaluate([10000.,3000.],method)
        original=SharedEntityLikelihood(['A']*3+['B']*3,{'signed':sparse.csr_matrix(z)},f,d,sx,sy).evaluate([10000.,3000.],method)
        np.testing.assert_allclose(spectral['negative_profiled_likelihood'],original['negative_profiled_likelihood'],rtol=2e-9,atol=2e-9)
        np.testing.assert_allclose(spectral['gradient'],original['gradient'],rtol=2e-6,atol=2e-10)
    # Use exactly the established independent-optimizer fixture, not production data.
    rng=np.random.default_rng(473);n=36;labels=np.repeat(np.arange(3),12);operators={}
    for name in ['gene','family']:
        z=np.zeros((n,9))
        for column in range(9):z[labels==column%3,column]=rng.normal(size=12)
        operators[name]=sparse.csr_matrix(z)
    f=rng.normal(size=(n,3));x=np.column_stack([np.ones(n),rng.normal(size=n)]);d=np.ones(n)
    y=operators['gene']@rng.normal(size=9)+.3*(operators['family']@rng.normal(size=9))+rng.normal(size=n)
    numerical=SharedEntityLikelihood(labels,operators,f,d,x,y)
    independent=IndependentSharedEntityLikelihood(labels,operators,f,d,x,y,column_batch=2)
    checked_candidates=0;curvature_checks=[];maximum_curvature_difference=0.
    for method in ['ml','reml']:
        candidate=fit_shared_entity(numerical,method,gradient_tolerance=3e-6)
        audit=audit_candidate(independent,candidate,gradient_atol=3e-6)
        assert audit['status']=='independent_shared_entity_replay_passed_pending_curvature_and_calibration',audit
        curvature=audit_local_curvature(independent,candidate,gradient_atol=3e-6)
        assert curvature['status']=='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits',curvature
        assert curvature['curvature_evaluations']==1+4*len(candidate['variance_ratios'])
        # Independent dense score derivatives at the finer bounded stencil.
        norms=np.asarray(candidate['kernel_normalization']);coordinates=np.log1p(np.asarray(candidate['variance_ratios'])*norms)
        upper=np.log1p(candidate['maximum_scaled_variance']);hessian=np.empty((len(norms),len(norms)))
        def dense_score(point):return np.asarray(dense(labels,operators,f,d,x,y,np.expm1(point)/norms,method)['gradient'])*np.exp(point)/norms
        initial=dense_score(coordinates)
        for j,value in enumerate(coordinates):
            h=min(5e-5*max(1.,value),upper/4)
            def change(distance):
                point=coordinates.copy();point[j]+=distance*h;return dense_score(point)
            if value>=h and upper-value>=h:hessian[:,j]=(change(1)-change(-1))/(2*h)
            elif upper-value>=2*h:hessian[:,j]=(-3*initial+4*change(1)-change(2))/(2*h)
            else:hessian[:,j]=(3*initial-4*change(-1)+change(-2))/(2*h)
        hessian=(hessian+hessian.T)/2
        np.testing.assert_allclose(curvature['score_hessian'],hessian,rtol=3e-5,atol=2e-7)
        maximum_curvature_difference=max(maximum_curvature_difference,float(np.max(abs(np.asarray(curvature['score_hessian'])-hessian))))
        curvature_checks.append({k:curvature[k] for k in ['status','curvature_evaluations','strong_lower_boundary_indices',
            'tested_curvature_eigenvalues','hessian_antisymmetry_norm','two_resolution_hessian_difference_norm','numerical_curvature_review_envelope']})
        checked_candidates+=1
    altered=[]
    def false_export(name,change):
        false=copy.deepcopy(candidate);change(false)
        try:result=audit_candidate(independent,false,gradient_atol=3e-6)
        except (ValueError,ArithmeticError):altered.append(name);return
        assert result['status']=='independent_shared_entity_replay_requires_review',(name,result)
        altered.append(name)
    false_export('objective',lambda c:c.update(negative_profiled_likelihood=c['negative_profiled_likelihood']+1.))
    false_export('nonfinite_objective',lambda c:c.update(negative_profiled_likelihood=float('nan')))
    false_export('coefficient',lambda c:c['beta'].__setitem__(0,c['beta'][0]+1.))
    false_export('conditional_covariance',lambda c:c['conditional_beta_covariance'][0].__setitem__(0,c['conditional_beta_covariance'][0][0]+1.))
    false_export('profiled_scale',lambda c:c.update(profiled_scale=c['profiled_scale']*2))
    false_export('variance_components',lambda c:c['variance_components'].__setitem__(0,c['variance_components'][0]+1.))
    false_export('scientific_acceptance_flag',lambda c:c.update(scientific_eligibility=True))
    false_export('conditional_covariance_flag',lambda c:c.update(coefficient_covariance_is_conditional=False))
    false_export('parameter_order',lambda c:c.update(parameter_names=list(reversed(c['parameter_names']))))
    false_export('kernel_normalization',lambda c:c['kernel_normalization'].__setitem__(0,c['kernel_normalization'][0]*2))
    false_export('negative_ratio',lambda c:c['variance_ratios'].__setitem__(0,-1.))
    false_export('unqualified_optimizer',lambda c:c.update(status='shared_entity_optimizer_candidate_requires_review'))
    false_export('false_variance_cap',lambda c:c.update(maximum_scaled_variance=1e-10))
    units=np.array([1e100,1e-100]);scaled_candidate=copy.deepcopy(candidate)
    scaled_candidate['beta']=(np.asarray(candidate['beta'])/units).tolist()
    scaled_candidate['conditional_beta_covariance']=(np.asarray(candidate['conditional_beta_covariance'])/units[:,None]/units[None,:]).tolist()
    scaled_independent=IndependentSharedEntityLikelihood(labels,operators,f,d,x*units,y)
    assert audit_candidate(scaled_independent,scaled_candidate,gradient_atol=3e-6)['status']=='independent_shared_entity_replay_passed_pending_curvature_and_calibration'
    for field in ['coefficient','conditional_covariance']:
        corrupted=copy.deepcopy(scaled_candidate)
        if field=='coefficient':corrupted['beta'][0]+=1e-98
        else:corrupted['conditional_beta_covariance'][0][0]+=1e-198
        assert audit_candidate(scaled_independent,corrupted,gradient_atol=3e-6)['status']=='independent_shared_entity_replay_requires_review'
        altered.append('tiny_original_unit_'+field)
    false=copy.deepcopy(candidate);false['negative_profiled_likelihood']+=1.
    failed_curvature=audit_local_curvature(independent,false)
    assert failed_curvature['status']=='independent_shared_entity_curvature_requires_review' and failed_curvature['curvature_evaluations']==0
    # Exact analytic score fields exercise constrained-curvature dispositions;
    # these are separate software stress fixtures, not covariance models.
    norms=np.asarray(candidate['kernel_normalization'])
    class AnalyticCurvatureFixture(IndependentSharedEntityLikelihood):
        def __init__(self,hessian,linear):
            super().__init__(labels,operators,f,d,x,y)
            self.hessian=hessian;self.linear=linear
        def evaluate(self,ratios,method='reml'):
            ratios=np.asarray(ratios);point=np.log1p(ratios*norms)
            return dict(negative_profiled_likelihood=1.+float(self.linear@point+.5*point@self.hessian@point),
                beta=np.zeros(2),conditional_beta_covariance=np.eye(2),profiled_scale=1.,inverse_relative_residual=0.,
                gradient=(self.linear+self.hessian@point)*norms/np.exp(point))
    analytic_candidate=copy.deepcopy(candidate)
    analytic_candidate.update(variance_ratios=[0.,0.,0.],variance_components=[0.,0.,0.],profiled_scale=1.,
        beta=[0.,0.],conditional_beta_covariance=np.eye(2).tolist(),negative_profiled_likelihood=1.)
    critical_cone_checks=[]
    for name,hessian,linear,expected in [
        ('weak_lower_positive',np.diag([2.,3.,4.]),np.zeros(3),True),
        ('weak_lower_negative',np.diag([-2.,3.,4.]),np.zeros(3),False),
        ('weak_lower_flat',np.zeros((3,3)),np.zeros(3),False),
        ('one_strong_lower',np.diag([-2.,3.,4.]),np.array([1.,0.,0.]),True),
        ('all_strong_lower',-np.eye(3),np.ones(3),True),
        ('inconsistent_asymmetric_score',np.array([[2.,1.,0.],[0.,3.,0.],[0.,0.,4.]]),np.zeros(3),False)]:
        result=audit_local_curvature(AnalyticCurvatureFixture(hessian,linear),analytic_candidate)
        assert (result['status']=='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits')==expected,(name,result)
        critical_cone_checks.append(dict(fixture=name,expected_pass=expected,status=result['status'],
            strong_lower_boundary_indices=result['strong_lower_boundary_indices']))
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,ArithmeticError):rejected.append(name);return
        raise AssertionError('Invalid independent input accepted: '+name)
    reject('cross_component_entity',lambda:IndependentSharedEntityLikelihood(np.arange(n),{'z':sparse.csr_matrix(np.ones((n,1)))},f,d,x,y))
    reject('rank_deficient_design',lambda:IndependentSharedEntityLikelihood(labels,operators,f,d,np.ones((n,2)),y))
    reject('perfect_response_fit',lambda:IndependentSharedEntityLikelihood(labels,operators,f,d,x,x@np.array([1.,2.])).evaluate([1.,1.,1.]))
    reject('negative_diagonal',lambda:IndependentSharedEntityLikelihood(labels,operators,f,-d,x,y))
    reject('reserved_species_name',lambda:IndependentSharedEntityLikelihood(labels,{'species':operators['gene']},f,d,x,y))
    reject('zero_column_batch',lambda:IndependentSharedEntityLikelihood(labels,operators,f,d,x,y,column_batch=0))
    reject('extreme_component_condition',lambda:IndependentSharedEntityLikelihood(labels,operators,f,d,x,y).evaluate([1e18,0.,0.]))
    result=dict(status='passed_independent_spectral_shared_entity_likelihood_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        independent_dense_ml_reml_cases=cases,maximum_absolute_objective_error=maximum_objective,
        maximum_absolute_gradient_error=maximum_gradient,maximum_inverse_relative_residual=maximum_inverse,
        strong_signed_cases=2,precision_digits=80,strong_signed_high_precision_errors=precision_errors,
        independent_optimizer_candidates_replayed=checked_candidates,altered_candidate_exports_rejected=altered,
        independent_local_curvature_checks=curvature_checks,maximum_dense_curvature_absolute_difference=maximum_curvature_difference,
        failed_candidate_curvature_refused=True,analytic_critical_cone_stress_fixtures=critical_cone_checks,
        invalid_inputs_rejected=rejected,entity_batch_bound_verified=2,permutations_and_extreme_design_units_passed=True,
        source_hashes={str(p.relative_to(Path.cwd()) if p.is_absolute() else p):sha(p) for p in [Path(__file__),Path('scripts/independent_shared_entity_likelihood.py'),
            Path('scripts/shared_entity_likelihood.py'),Path('scripts/shared_entity_covariance.py'),
            Path('scripts/fit_shared_entity_likelihood.py'),Path('scripts/check_shared_entity_likelihood.py'),
            Path('scripts/covariance_basis_context.py'),Path('scripts/covariance_basis_audit.py')]},
        scope='Independent component eigendecomposition, species SVD and fixed-design SVD; '
            'streamed positive-norm kernel traces. Synthetic software fixtures only. No production '
            'candidate acceptance, rigorous curvature error bound/global optimum or calibrated inference. '
            'Two covariance-fit curvature replays and six explicitly analytic critical-cone stress fixtures.')
    with a.output.open('x') as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
