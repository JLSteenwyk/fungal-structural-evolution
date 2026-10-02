"""Profile ML/REML with analytic variance-ratio scores for explicit kernels.

Numerical fitting support only. Covariance definitions, cohort identifiability,
adequacy, multiple testing and calibration must be qualified separately.
The launched shared-entity operator implementation remains unchanged.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve,qr,solve_triangular
from shared_entity_covariance import SharedEntityCovariance


class SharedEntityLikelihood:
    def __init__(self,labels,incidence,factor,diagonal,design,response):
        self.labels=np.asarray(labels);self.incidence=incidence
        self.factor=np.asarray(factor,dtype=float);self.diagonal=np.asarray(diagonal,dtype=float)
        self.y=np.asarray(response,dtype=float);x=np.asarray(design,dtype=float);self.design=x.copy()
        n=len(self.labels)
        if x.ndim!=2 or x.shape[0]!=n or self.y.shape!=(n,) or not 0<x.shape[1]<n:
            raise ValueError('Compatible active design, response and residual dimensions required')
        if not np.isfinite(x).all() or not np.isfinite(self.y).all():raise ValueError('Finite design and response required')
        self.maximum=np.max(abs(x),axis=0)
        if np.any(self.maximum==0):raise ValueError('Inactive zero design columns require explicit separate accounting')
        scaled=x/self.maximum;self.length=np.sqrt(np.sum(scaled*scaled,axis=0));scaled/=self.length
        singular=np.linalg.svd(scaled,compute_uv=False)
        if singular[-1]<=10*max(x.shape)*np.finfo(float).eps*singular[0]:raise ValueError('Fixed design requires rank/boundary review')
        self.q,self.r=qr(scaled,mode='economic');self.n=n;self.p=x.shape[1]
        self.design_log_jacobian=float(2*(np.log(self.maximum).sum()+np.log(self.length).sum()+np.log(abs(np.diag(self.r))).sum()))
        # Use the existing strict covariance validation, even for variance-zero
        # columns, and retain its exact component and coalescence definitions.
        template=SharedEntityCovariance(self.labels,incidence,self.factor,self.diagonal,{k:0. for k in incidence},0.)
        self.incidence=template.incidence;self.parts=template.blocks
        self.kernels=[{name:(z[rows]@z[rows].T).toarray() for name,z in self.incidence.items()} for rows in self.parts]
        self.parameter_names=[*self.incidence,'species'];self.last_covariance=None

    def evaluate(self,ratios,method='reml',gradient=True):
        ratios=np.asarray(ratios,dtype=float)
        if ratios.shape!=(len(self.parameter_names),) or not np.isfinite(ratios).all() or np.any(ratios<0):
            raise ValueError('One finite nonnegative variance ratio per named kernel required')
        if method not in ['ml','reml']:raise ValueError('Likelihood must be ml or reml')
        variances=dict(zip(self.incidence,ratios[:-1]));cov=SharedEntityCovariance(self.labels,self.incidence,self.factor,self.diagonal,variances,ratios[-1])
        solved=cov.solve(np.column_stack([self.q,self.y]));iq=solved[:,:self.p];iy=solved[:,-1]
        info=self.q.T@iq;info=(info+info.T)/2;lower=cho_factor(info,lower=True,check_finite=True)
        coefficient=cho_solve(lower,self.q.T@iy,check_finite=False)
        residual=self.y-self.q@coefficient;alpha=cov.solve(residual)
        quadratic=float(residual@alpha);reference=float(self.y@iy)
        if not np.isfinite(quadratic) or quadratic<=64*np.finfo(float).eps**2*max(reference,np.finfo(float).tiny):
            raise ArithmeticError('Residual quadratic is zero or below numerical resolution')
        df=self.n-self.p if method=='reml' else self.n;scale=quadratic/df
        determinant=float(cov.logdet)
        if method=='reml':determinant+=float(2*np.log(np.diag(lower[0])).sum())+self.design_log_jacobian
        objective=.5*(determinant+df*(1+np.log(2*np.pi*scale)))
        inverse_info=cho_solve(lower,np.eye(self.p),check_finite=False)
        inverse_r=solve_triangular(self.r,np.eye(self.p),lower=False,check_finite=False)
        beta=solve_triangular(self.r,coefficient,lower=False,check_finite=False)/self.length/self.maximum
        beta_covariance=(inverse_r@inverse_info@inverse_r.T)*scale
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            beta_covariance=beta_covariance/self.length[:,None]/self.maximum[:,None]/self.length[None,:]/self.maximum[None,:]
        components=scale*ratios
        if not np.isfinite(objective) or not np.isfinite(beta).all() or not np.isfinite(beta_covariance).all() or np.any(np.diag(beta_covariance)<=0) or not np.isfinite(components).all():
            raise ArithmeticError('Likelihood or original-unit coefficients are not numerically representable')
        result=dict(method=method,negative_profiled_likelihood=float(objective),beta=beta,
            conditional_beta_covariance=beta_covariance,profiled_scale=float(scale),residual_quadratic=quadratic,
            scale_degrees_of_freedom=df,parameter_names=self.parameter_names,
            variance_ratios=ratios.copy(),variance_components=components,
            coefficient_covariance_is_conditional=True,likelihood_evaluations=1)
        self.last_covariance=cov
        if not gradient:return result
        traces={name:0. for name in self.incidence}
        for rows,cholesky,kernels in zip(self.parts,cov.cholesky,self.kernels):
            inverse=cho_solve(cholesky,np.eye(len(rows)),check_finite=False)
            if cov.phylogenetic_cholesky is not None:
                projection=cov.base_inverse_factor[rows]
                inverse-=cov.pvar*(projection@cho_solve(cov.phylogenetic_cholesky,projection.T,check_finite=False))
            if not np.isfinite(inverse).all() or np.any(np.diag(inverse)<=0):
                raise ArithmeticError('Conditional block inverse lost finite positivity')
            for name,kernel in kernels.items():traces[name]+=float(np.sum(inverse*kernel))
        species_inverse=cov.solve(self.factor)
        species_trace=float(np.sum(self.factor*species_inverse))
        scores=[];trace_values=[];energy_values=[]
        for name in self.parameter_names:
            if name=='species':
                trace=species_trace;projected=self.factor.T@iq;energy=float(np.sum((self.factor.T@alpha)**2))
            else:
                z=self.incidence[name];trace=traces[name];projected=z.T@iq;energy=float(np.sum((z.T@alpha)**2))
            if method=='reml':trace-=float(np.sum(inverse_info*(projected.T@projected)))
            scores.append(.5*(trace-df*energy/quadratic));trace_values.append(trace);energy_values.append(energy)
        if not np.isfinite(scores).all():raise ArithmeticError('Nonfinite analytic likelihood gradient')
        result.update(gradient=np.asarray(scores),kernel_traces_in_score=np.asarray(trace_values),
            residual_kernel_energies=np.asarray(energy_values))
        return result
