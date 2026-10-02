#!/usr/bin/env python3
"""Dense and finite-difference checks of ML/REML analytic variance scores."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from shared_entity_likelihood import SharedEntityLikelihood
from run_ortholog_pair_guide_comparison import sha


def dense(labels,incidence,factor,diagonal,x,y,ratios,method):
    kernels=[z.toarray()@z.toarray().T for z in incidence.values()]+[factor@factor.T]
    v=np.diag(diagonal)+sum(r*k for r,k in zip(ratios,kernels));iv=np.linalg.inv(v)
    information=x.T@iv@x;ii=np.linalg.inv(information);beta=ii@x.T@iv@y
    residual=y-x@beta;alpha=iv@residual;quadratic=residual@alpha
    projection=iv-iv@x@ii@x.T@iv;n,p=x.shape;df=n-p if method=='reml' else n
    objective=.5*(np.linalg.slogdet(v)[1]+(np.linalg.slogdet(information)[1] if method=='reml' else 0)
        +df*(1+np.log(2*np.pi*quadratic/df)))
    gradients=[.5*(np.trace((projection if method=='reml' else iv)@k)-df*(alpha@k@alpha)/quadratic) for k in kernels]
    return dict(negative_profiled_likelihood=objective,beta=beta,conditional_beta_covariance=ii*(quadratic/df),
        profiled_scale=quadratic/df,gradient=gradients,residual_quadratic=quadratic)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(328);cases=0;max_objective=max_gradient=max_numeric=0.
    for trial in range(18):
        n=12+trial%9;labels=np.arange(n)%3;operators={}
        for name in ['gene','model','family']:
            z=np.zeros((n,9))
            for column in range(9):z[labels==column%3,column]=rng.integers(-2,3,size=(labels==column%3).sum())/2
            operators[name]=sparse.csr_matrix(z)
        if trial%3==0:operators['cancelled_family']=sparse.csr_matrix((n,5))
        f=rng.normal(size=(n,trial%5));d=rng.uniform(.3,2.,n)
        x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n)
        ratios=rng.uniform(.01,2.,len(operators)+1)
        if trial%2==0:ratios[0]=0.
        if trial%6==0:ratios[:]=0.
        likelihood=SharedEntityLikelihood(labels,operators,f,d,x,y)
        for method in ['ml','reml']:
            actual=likelihood.evaluate(ratios,method);expected=dense(labels,operators,f,d,x,y,ratios,method)
            for key in expected:np.testing.assert_allclose(actual[key],expected[key],rtol=3e-9,atol=3e-10)
            max_objective=max(max_objective,abs(actual['negative_profiled_likelihood']-expected['negative_profiled_likelihood']))
            max_gradient=max(max_gradient,float(np.max(abs(actual['gradient']-expected['gradient']))))
            for j,value in enumerate(ratios):
                h=1e-4*max(1.,value)
                def point(step):
                    changed=ratios.copy();changed[j]+=step*h
                    return likelihood.evaluate(changed,method,False)['negative_profiled_likelihood']
                numeric=((-25*point(0)+48*point(1)-36*point(2)+16*point(3)-3*point(4))/(12*h)
                    if value<2*h else (point(-2)-8*point(-1)+8*point(1)-point(2))/(12*h))
                np.testing.assert_allclose(actual['gradient'][j],numeric,rtol=3e-5,atol=3e-7)
                max_numeric=max(max_numeric,abs(actual['gradient'][j]-numeric))
            permutation=rng.permutation(n)
            permuted=SharedEntityLikelihood(labels[permutation],{k:z[permutation] for k,z in operators.items()},f[permutation],d[permutation],x[permutation],y[permutation]).evaluate(ratios,method)
            for key in expected:np.testing.assert_allclose(actual[key],permuted[key],rtol=3e-9,atol=3e-10)
            cases+=1
    # The same fixed column space in very different units must preserve scores.
    n=16;labels=np.zeros(n);z={'shared':sparse.csr_matrix(rng.normal(size=(n,4)))};f=rng.normal(size=(n,2))
    x=np.column_stack([np.ones(n),rng.normal(size=n)]);y=rng.normal(size=n);ratios=np.array([1.2,.3]);units=np.array([1e100,1e-100])
    base=SharedEntityLikelihood(labels,z,f,np.ones(n),x,y).evaluate(ratios)
    scaled=SharedEntityLikelihood(labels,z,f,np.ones(n),x*units,y).evaluate(ratios)
    np.testing.assert_allclose(base['gradient'],scaled['gradient'],rtol=3e-9,atol=3e-10)
    np.testing.assert_allclose(base['beta'],scaled['beta']*units,rtol=3e-9,atol=3e-10)
    np.testing.assert_allclose(base['conditional_beta_covariance'],scaled['conditional_beta_covariance']*units[:,None]*units[None,:],rtol=3e-9,atol=3e-10)
    invalid=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,ArithmeticError):invalid.append(name)
        else:raise AssertionError('Invalid likelihood input accepted: '+name)
    reject('negative_ratio',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,y).evaluate([-1,.3]))
    reject('nonfinite_ratio',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,y).evaluate([float('nan'),.3]))
    reject('wrong_ratio_count',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,y).evaluate([1.]))
    reject('unknown_method',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,y).evaluate(ratios,'other'))
    reject('rank_deficient_design',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),np.ones((n,2)),y))
    reject('nonfinite_response',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,np.full(n,np.inf)))
    reject('perfect_response_fit',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x,x@np.array([1.,2.])).evaluate(ratios))
    reject('unrepresentable_coefficient_covariance',lambda:SharedEntityLikelihood(labels,z,f,np.ones(n),x*np.array([1e-200,1.]),y).evaluate(ratios))
    reject('cross_component_entity',lambda:SharedEntityLikelihood(np.arange(n),{'z':sparse.csr_matrix(np.ones((n,1)))},f,np.ones(n),x,y))
    result=dict(status='passed_shared_entity_analytic_profile_likelihood_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        ml_reml_dense_cases=cases,maximum_absolute_objective_error=max_objective,maximum_absolute_analytic_gradient_error=max_gradient,
        maximum_absolute_finite_difference_error=max_numeric,row_permutation_passed=True,extreme_design_units_passed=True,
        zero_variance_boundaries_and_zero_rank_species_passed=True,rejected_invalid_inputs=invalid,
        source_hashes={str(p):sha(p) for p in [Path(__file__),Path('scripts/shared_entity_likelihood.py'),Path('scripts/shared_entity_covariance.py')]},
        scope='Numerical likelihood/gradient contracts only; dense inverses and forward/central fourth-order '
              'numerical derivatives. Not an optimizer, production fit, a biological pilot or uncertainty calibration.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
