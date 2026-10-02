#!/usr/bin/env python3
"""Compare shared-entity operators to independent dense and high-precision algebra."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy import sparse
from shared_entity_covariance import SharedEntityCovariance, profiled_reml
from run_ortholog_pair_guide_comparison import sha


def dense_covariance(diagonal,operators,factor,variances,species):
    n=len(diagonal);v=np.diag(diagonal)
    for name,z in operators.items():
        # Explicit sum of each latent entity's outer product, not sparse Gram.
        for j in range(z.shape[1]):
            u=np.asarray(z[:,j].toarray()).ravel();v+=variances[name]*np.outer(u,u)
    for j in range(factor.shape[1]):v+=species*np.outer(factor[:,j],factor[:,j])
    return v


def dense_reml(v,x,y):
    i=np.linalg.inv(v);information=x.T@i@x;inverse=np.linalg.inv(information);beta=inverse@x.T@i@y
    residual=y-x@beta;q=float(residual@i@residual);n,p=x.shape;scale=q/(n-p)
    return dict(negative_profiled_reml=.5*(np.linalg.slogdet(v)[1]+np.linalg.slogdet(information)[1]+(n-p)*(1+np.log(2*np.pi*scale))),
        beta=beta,profiled_scale=scale,residual_quadratic=q,conditional_beta_covariance=scale*inverse,residual_degrees_of_freedom=n-p)


def rejected(function):
    try:function()
    except (ValueError,ArithmeticError):return
    raise AssertionError('Invalid covariance/design accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rng=np.random.default_rng(309);checked=0;maximum_solve=maximum_logdet=maximum_objective=0.
    for trial in range(36):
        n=8+trial%17;blocks=np.arange(n)%3;factor=rng.normal(size=(n,trial%7));operators={}
        for name in ['pair','gene','model','family_intercept']:
            z=np.zeros((n,9))
            for j in range(9):z[blocks==j%3,j]=rng.integers(-2,3,size=int((blocks==j%3).sum()))/2
            operators[name]=sparse.csr_matrix(z)
        # Signed family endpoints cancel; preserve an explicit zero basis.
        operators['signed_family']=sparse.csr_matrix((n,3))
        d=rng.uniform(.5,2,size=n);variances={name:(0. if trial%3==0 else rng.uniform(0,3)) for name in operators}
        sp=0. if trial%4==0 else .7
        op=SharedEntityCovariance(blocks,operators,factor,d,variances,sp);v=dense_covariance(d,operators,factor,variances,sp)
        rhs=rng.normal(size=(n,4));solution=op.solve(rhs);expected=np.linalg.solve(v,rhs)
        np.testing.assert_allclose(solution,expected,rtol=2e-10,atol=2e-11)
        np.testing.assert_allclose(op.solve(rhs[:,0]),expected[:,0],rtol=2e-10,atol=2e-11)
        np.testing.assert_allclose(op.apply(rhs),v@rhs,rtol=2e-12,atol=2e-12)
        np.testing.assert_allclose(op.inverse_diagonal(3),np.diag(np.linalg.inv(v)),rtol=2e-10,atol=2e-11)
        np.testing.assert_allclose(op.logdet,np.linalg.slogdet(v)[1],rtol=2e-11,atol=2e-11)
        x=np.column_stack([np.ones(n),rng.normal(size=(n,2))]);y=rng.normal(size=n)
        result=profiled_reml(op,x,y);reference=dense_reml(v,x,y)
        for name in reference:np.testing.assert_allclose(result[name],reference[name],rtol=3e-10,atol=3e-10)
        perm=rng.permutation(n);other=SharedEntityCovariance(blocks[perm],{k:z[perm] for k,z in operators.items()},factor[perm],d[perm],variances,sp)
        np.testing.assert_allclose(other.solve(rhs[perm]),solution[perm],rtol=3e-10,atol=3e-11)
        np.testing.assert_allclose(other.logdet,op.logdet,rtol=2e-11,atol=2e-11)
        maximum_solve=max(maximum_solve,float(np.max(abs(solution-expected))))
        maximum_logdet=max(maximum_logdet,abs(op.logdet-np.linalg.slogdet(v)[1]))
        maximum_objective=max(maximum_objective,abs(result['negative_profiled_reml']-reference['negative_profiled_reml']));checked+=1
    # Duplicated COO coordinates and opposite signs must coalesce, never create noise.
    dup=sparse.coo_matrix(([.5,-.5,1.,1.],([0,0,1,2],[0,0,1,1])),shape=(4,2))
    aliases=SharedEntityCovariance(['A']*4,{'model':dup},np.zeros((4,0)),np.ones(4),{'model':2.})
    np.testing.assert_array_equal(aliases.incidence['model'].toarray(),[[0,0],[0,1],[0,1],[0,0]])
    np.testing.assert_allclose(aliases.solve(np.ones(4)),[1,.2,.2,1],atol=1e-15)
    # Independent 80-digit inverse/determinant for strongly correlated signed inputs.
    d=np.array([.2,.5,.7,.3,.6,.4]);z=sparse.csr_matrix([[1,0],[1,0],[-1,0],[0,1],[0,1],[0,-1]],dtype=float)
    f=np.column_stack([np.linspace(.1,.6,6),np.array([.2,-.1,.3,-.2,.5,.4])]);ratios={'shared':1e4}
    op=SharedEntityCovariance(['A']*3+['B']*3,{'shared':z},f,d,ratios,3000.)
    v=dense_covariance(d,{'shared':z},f,ratios,3000.);rhs=np.array([.3,-.8,.1,.6,-.3,.2])
    with mp.workdps(80):
        md=mp.diag([mp.mpf(str(v)) for v in d]);mz=mp.matrix(z.toarray().tolist());mf=mp.matrix([[mp.mpf(str(v)) for v in row] for row in f])
        mv=md+mp.mpf('10000')*mz*mz.T+mp.mpf('3000')*mf*mf.T
        target=np.array(list(mp.lu_solve(mv,mp.matrix([mp.mpf(str(v)) for v in rhs]))),dtype=float).ravel()
        ld=float(mp.log(mp.det(mv)))
    np.testing.assert_allclose(op.solve(rhs),target,rtol=2e-8,atol=2e-10)
    np.testing.assert_allclose(op.logdet,ld,rtol=2e-10,atol=2e-10)
    # Invalid partitioning is rejected even when its variance currently equals zero.
    invalid=[]
    def test(name,fn):rejected(fn);invalid.append(name)
    eye=sparse.eye(4,format='csr');factor=np.zeros((4,0));labels=['A']*4
    test('entity_crosses_blocks',lambda:SharedEntityCovariance(['A','B','A','B'],{'z':sparse.csr_matrix(np.ones((4,1)))},factor,np.ones(4),{'z':0.}))
    for name,diag,var,sp in [('nonpositive_diagonal',[1,1,0,1],1,0),('nonfinite_diagonal',[1,1,float('nan'),1],1,0),('negative_entity',np.ones(4),-1,0),('nonfinite_entity',np.ones(4),float('inf'),0),('negative_species',np.ones(4),1,-1)]:
        test(name,lambda diag=diag,var=var,sp=sp:SharedEntityCovariance(labels,{'z':eye},factor,diag,{'z':var},sp))
    test('missing_named_variance',lambda:SharedEntityCovariance(labels,{'z':eye},factor,np.ones(4),{}))
    test('nonfinite_factor',lambda:SharedEntityCovariance(labels,{'z':eye},np.full((4,1),np.nan),np.ones(4),{'z':1.}))
    test('nonfinite_incidence',lambda:SharedEntityCovariance(labels,{'z':sparse.csr_matrix(np.full((4,1),np.nan))},factor,np.ones(4),{'z':1.}))
    test('nonfinite_rhs',lambda:aliases.solve(np.full(4,np.nan)))
    test('rank_deficient_design',lambda:profiled_reml(aliases,np.ones((4,2)),np.arange(4)))
    test('exact_zero_response',lambda:profiled_reml(aliases,np.ones((4,1)),np.zeros(4)))
    test('exact_nonzero_fit',lambda:profiled_reml(aliases,np.ones((4,1)),np.ones(4)))
    result=dict(status='passed_shared_entity_block_low_rank_covariance_contracts',dense_cases=checked,
        maximum_absolute_solve_error=maximum_solve,maximum_absolute_logdet_error=maximum_logdet,maximum_absolute_profiled_reml_error=maximum_objective,
        high_precision_signed_shared_case_passed=True,row_permutation_passed=True,zero_rank_and_zero_variance_passed=True,
        duplicate_and_cancelling_aliases_preserved=True,invalid_cases_rejected=invalid,
        source_hashes={str(p):sha(p) for p in [Path(__file__),Path('scripts/shared_entity_covariance.py')]},
        scope='Numerical software contracts only. Independent explicit latent-entity outer products/dense inverses and 80-digit inverse/determinant; no biological pilot, production fit, optimizer or inferential calibration.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
