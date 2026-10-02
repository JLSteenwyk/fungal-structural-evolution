#!/usr/bin/env python3
"""Exact block shared-entity covariance with a low-rank phylogenetic update.

Numerical operator only: supplied loadings and variances are explicit working
model choices. Optimizer, adequacy, identifiability and calibration are separate.
"""
import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve


class SharedEntityCovariance:
    def __init__(self, block_labels, incidence, species_factor, residual_diagonal,
                 variances, species_variance=0.):
        labels=np.asarray(block_labels)
        if labels.ndim!=1 or not len(labels):raise ValueError('Nonempty one-dimensional block labels required')
        self.n=len(labels);_,self.block=np.unique(labels,return_inverse=True)
        self.diagonal=np.asarray(residual_diagonal,dtype=float)
        if self.diagonal.shape!=(self.n,) or not np.isfinite(self.diagonal).all() or np.any(self.diagonal<=0):
            raise ValueError('Residual diagonal must be positive finite with n rows')
        if set(variances)!=set(incidence):raise ValueError('Exactly one named variance per incidence operator required')
        self.variances={k:float(v) for k,v in variances.items()}
        if any(not np.isfinite(v) or v<0 for v in self.variances.values()):raise ValueError('Finite nonnegative entity variances required')
        self.pvar=float(species_variance)
        if not np.isfinite(self.pvar) or self.pvar<0:raise ValueError('Finite nonnegative species variance required')
        self.factor=np.asarray(species_factor,dtype=float)
        if self.factor.ndim!=2 or self.factor.shape[0]!=self.n or not np.isfinite(self.factor).all():
            raise ValueError('Finite n-by-r species factor required')
        self.incidence={}
        for name,original in incidence.items():
            z=sparse.csr_matrix(original,dtype=float,copy=True)
            if z.shape[0]!=self.n or not np.isfinite(z.data).all():raise ValueError('Invalid entity incidence operator')
            z.sum_duplicates();z.eliminate_zeros();z.sort_indices()
            rows,cols=z.nonzero();small=np.full(z.shape[1],len(labels),dtype=np.int64);large=np.full(z.shape[1],-1,dtype=np.int64)
            np.minimum.at(small,cols,self.block[rows]);np.maximum.at(large,cols,self.block[rows])
            present=large>=0
            if np.any(small[present]!=large[present]):raise ValueError('An entity spans multiple declared blocks')
            self.incidence[name]=z
        order=np.argsort(self.block,kind='stable');cuts=np.flatnonzero(np.diff(self.block[order]))+1
        self.blocks=np.split(order,cuts);self.cholesky=[];self.logdet=0.
        for rows in self.blocks:
            core=np.diag(self.diagonal[rows])
            for name,z in self.incidence.items():
                v=self.variances[name]
                if v:
                    local=z[rows];core+=v*(local@local.T).toarray()
            lower=cho_factor(core,lower=True,check_finite=True)
            self.cholesky.append(lower);self.logdet+=float(2*np.log(np.diag(lower[0])).sum())
        self.phylogenetic_cholesky=None
        if self.pvar>0 and self.factor.shape[1]:
            self.base_inverse_factor=self.base_solve(self.factor)
            core=np.eye(self.factor.shape[1])+self.pvar*(self.factor.T@self.base_inverse_factor)
            core=(core+core.T)/2
            self.phylogenetic_cholesky=cho_factor(core,lower=True,check_finite=True)
            self.logdet+=float(2*np.log(np.diag(self.phylogenetic_cholesky[0])).sum())

    def _values(self,values):
        v=np.asarray(values,dtype=float);vector=v.ndim==1
        if vector:v=v[:,None]
        if v.ndim!=2 or v.shape[0]!=self.n or not np.isfinite(v).all():raise ValueError('Finite right-hand side with n rows required')
        return v,vector

    def base_solve(self,values):
        values,vector=self._values(values);out=np.empty_like(values)
        for rows,lower in zip(self.blocks,self.cholesky):out[rows]=cho_solve(lower,values[rows],check_finite=False)
        return out[:,0] if vector else out

    def solve(self,values):
        values,vector=self._values(values);out=self.base_solve(values)
        if self.phylogenetic_cholesky is not None:
            correction=cho_solve(self.phylogenetic_cholesky,self.factor.T@out,check_finite=False)
            out-=self.pvar*(self.base_inverse_factor@correction)
        # Correct low-rank subtraction roundoff without changing the covariance.
        limit=1e-9*max(float(np.max(abs(values),initial=0)),np.finfo(float).tiny)
        for step in range(4):
            residual=values-self.apply(out);error=float(np.max(abs(residual),initial=0))
            if error<=limit:
                self.last_solve_diagnostic=dict(refinement_steps=step,maximum_residual=error,relative_residual_limit=limit)
                break
            if step==3:raise ArithmeticError('Covariance solve did not meet residual check after three refinements')
            update=self.base_solve(residual)
            if self.phylogenetic_cholesky is not None:
                update-=self.pvar*(self.base_inverse_factor@cho_solve(self.phylogenetic_cholesky,self.factor.T@update,check_finite=False))
            out+=update
        return out[:,0] if vector else out

    def apply(self,values):
        values,vector=self._values(values);out=self.diagonal[:,None]*values
        for name,z in self.incidence.items():
            if self.variances[name]:out+=self.variances[name]*(z@(z.T@values))
        if self.pvar:out+=self.pvar*(self.factor@(self.factor.T@values))
        return out[:,0] if vector else out

    def inverse_diagonal(self,block_columns=32):
        if not isinstance(block_columns,int) or block_columns<1:raise ValueError('Positive integer diagonal block size required')
        out=np.empty(self.n)
        for rows,lower in zip(self.blocks,self.cholesky):
            for start in range(0,len(rows),block_columns):
                stop=min(start+block_columns,len(rows));e=np.zeros((len(rows),stop-start));e[np.arange(start,stop),np.arange(stop-start)]=1
                solved=cho_solve(lower,e,check_finite=False);out[rows[start:stop]]=solved[np.arange(start,stop),np.arange(stop-start)]
        if self.phylogenetic_cholesky is not None:
            correction=cho_solve(self.phylogenetic_cholesky,self.base_inverse_factor.T,check_finite=False)
            out-=self.pvar*np.sum(self.base_inverse_factor*correction.T,axis=1)
        if not np.isfinite(out).all() or np.any(out<=0):raise ArithmeticError('Inverse diagonal lost numerical positivity')
        return out


def profiled_reml(covariance, design, response):
    """Direct-residual profiled REML and conditional coefficients, never a fit."""
    x=np.asarray(design,dtype=float);y=np.asarray(response,dtype=float)
    if x.ndim!=2 or x.shape[0]!=covariance.n or y.shape!=(covariance.n,):raise ValueError('Incompatible design/response')
    n,p=x.shape
    if not p or n<=p or not np.isfinite(x).all() or not np.isfinite(y).all() or np.linalg.matrix_rank(x)!=p:
        raise ValueError('Finite full-rank design and positive residual dimensions required')
    ix=covariance.solve(x);iy=covariance.solve(y);information=x.T@ix;information=(information+information.T)/2
    lower=cho_factor(information,lower=True,check_finite=True);beta=cho_solve(lower,x.T@iy)
    residual=y-x@beta;ir=covariance.solve(residual);q=float(residual@ir)
    reference=float(y@iy)
    if not np.isfinite(q) or q<=64*np.finfo(float).eps**2*max(reference,np.finfo(float).tiny):
        raise ArithmeticError('Residual quadratic form is zero or below its numerical resolution')
    df=n-p;scale=q/df;logdet_info=float(2*np.log(np.diag(lower[0])).sum())
    value=.5*(covariance.logdet+logdet_info+df*(1+np.log(2*np.pi*scale)))
    if not np.isfinite(value):raise ArithmeticError('Nonfinite profiled REML')
    return dict(negative_profiled_reml=float(value),beta=beta,profiled_scale=float(scale),residual_quadratic=q,
        conditional_beta_covariance=scale*cho_solve(lower,np.eye(p)),residual_degrees_of_freedom=df)
