"""Component-local spectral traces with explicit-whitening precision fallback.

Independent of production Cholesky/variance-fit arithmetic. Inherits the
validated input/boundary protocol; the original streamed checker stays intact.
"""
import numpy as np
from scipy import linalg
from independent_shared_entity_likelihood import IndependentSharedEntityLikelihood


class ComponentSpectralLikelihood(IndependentSharedEntityLikelihood):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.local_operators=[];self.component_kernels=[];widths=np.zeros(len(self.incidence),dtype=int)
        for rows in self.parts:
            operators=[];kernels=[]
            for j,z in enumerate(self.incidence.values()):
                local=z[rows];columns=np.unique(local.indices);local=local[:,columns].tocsr()
                operators.append((columns,local));widths[j]+=len(columns)
                kernels.append((local@local.T).toarray() if len(columns) else None)
            self.local_operators.append(operators);self.component_kernels.append(kernels)
        assert widths.tolist()==[z.shape[1] for z in self.incidence.values()]

    def evaluate(self,ratios,method='reml'):
        ratios=np.asarray(ratios,dtype=float)
        if ratios.shape!=(len(self.parameter_names),) or not np.isfinite(ratios).all() or np.any(ratios<0):
            raise ValueError('Explicit finite nonnegative ratios required')
        if method not in ['ml','reml']:raise ValueError('Likelihood must be ml or reml')
        roots=[];logdet=0.;maximum_condition=1.
        for rows,kernels in zip(self.parts,self.component_kernels):
            base=np.diag(self.diagonal[rows])
            for ratio,kernel in zip(ratios[:-1],kernels):
                if ratio and kernel is not None:base+=ratio*kernel
            if not np.isfinite(base).all():raise ArithmeticError('Unrepresentable component covariance')
            values,vectors=linalg.eigh(base,driver='evd',check_finite=True)
            if values[0]<=64*np.finfo(float).eps*len(rows)*values[-1]:
                raise ArithmeticError('Component inverse square root requires precision review')
            roots.append((vectors,1/np.sqrt(values)));logdet+=float(np.log(values).sum())
            maximum_condition=max(maximum_condition,float(values[-1]/values[0]))
        def local_root(values,root):
            vectors,inverse_root=root
            return vectors@(inverse_root[:,None]*(vectors.T@values))
        def base_root(values):
            out=np.empty_like(values)
            for rows,root in zip(self.parts,roots):out[rows]=local_root(values[rows],root)
            return out
        if ratios[-1] and self.factor.shape[1]:
            f=base_root(self.factor)*np.sqrt(ratios[-1])
            if not np.isfinite(f).all():raise ArithmeticError('Unrepresentable whitened species factor')
            u,s,_=linalg.svd(f,full_matrices=False,lapack_driver='gesvd',check_finite=True)
            species_root_values=1/np.hypot(1.,s);logdet+=float(2*np.log(np.hypot(1.,s)).sum())
        else:u=np.empty((self.n,0));species_root_values=np.empty(0)
        def species_root(values):return values+u@((species_root_values-1)[:,None]*(u.T@values))
        def whiten(values):return species_root(base_root(values))
        def inverse(values):return base_root(species_root(whiten(values)))
        def apply(values):
            out=self.diagonal[:,None]*values
            for ratio,z in zip(ratios[:-1],self.incidence.values()):
                if ratio:out+=ratio*(z@(z.T@values))
            if ratios[-1]:out+=ratios[-1]*(self.factor@(self.factor.T@values))
            return out
        rhs=np.column_stack([self.scaled,self.response]);solved=inverse(rhs)
        relative_error=float(np.max(abs(apply(solved)-rhs),initial=0))/max(float(np.max(abs(rhs),initial=0)),np.finfo(float).tiny)
        if not np.isfinite(relative_error) or relative_error>1e-8:raise ArithmeticError('Independent spectral inverse failed residual check')
        wx=whiten(self.scaled);wy=whiten(self.response[:,None])[:,0]
        left,values,right=linalg.svd(wx,full_matrices=False,lapack_driver='gesvd',check_finite=True)
        if values[-1]<=10*max(wx.shape)*np.finfo(float).eps*values[0]:raise ArithmeticError('Whitened fixed design requires precision review')
        coefficients=right.T@((left.T@wy)/values);residual=wy-wx@coefficients;quadratic=float(residual@residual)
        if not np.isfinite(quadratic) or quadratic<=64*np.finfo(float).eps**2*max(float(wy@wy),np.finfo(float).tiny):
            raise ArithmeticError('Residual quadratic is zero or below numerical resolution')
        df=self.n-self.p if method=='reml' else self.n;scale=quadratic/df
        if method=='reml':logdet+=float(2*np.log(values).sum())+self.column_log_jacobian
        objective=.5*(logdet+df*(1+np.log(2*np.pi*scale)));beta=coefficients/self.length/self.maximum
        covariance=(right.T/(values*values))@right*scale
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            covariance=covariance/self.length[:,None]/self.maximum[:,None]/self.length[None,:]/self.maximum[None,:]
        if not np.isfinite(objective) or not np.isfinite(beta).all() or not np.isfinite(covariance).all() or np.any(np.diag(covariance)<=0):
            raise ArithmeticError('Original-unit coefficients are not representable')
        fixed_species=(left.T@u)*(species_root_values-1)
        residual_species=(species_root_values-1)*(u.T@residual)
        species_penalty_weights=1-species_root_values**2
        count=len(self.incidence);traces=np.zeros(count+1);energies=np.zeros(count+1)
        local_batches=fallbacks=global_entity_cells=local_entity_cells=0;largest_batch=0
        # Each entity column belongs to exactly one verified component. All
        # contractions except a precision fallback stay within its row block.
        for rows,root,operators in zip(self.parts,roots,self.local_operators):
            local_u=u[rows];local_left=left[rows];local_residual=residual[rows]
            for j,(column_ids,z) in enumerate(operators):
                for start in range(0,z.shape[1],self.column_batch):
                    columns=z[:,start:start+self.column_batch].toarray();a=local_root(columns,root)
                    projected=local_u.T@a
                    base_norm=float(np.sum(a*a));species_penalty=float(np.sum(species_penalty_weights[:,None]*(projected*projected)))
                    fixed_image=local_left.T@a+fixed_species@projected
                    fixed_penalty=float(np.sum(fixed_image*fixed_image)) if method=='reml' else 0.
                    trace=base_norm-species_penalty-fixed_penalty
                    local_energy=a.T@local_residual;correction=projected.T@residual_species
                    energy_image=local_energy+correction
                    dimension=max(self.n,self.p,u.shape[1],len(rows),self.column_batch,1)
                    envelope=64*np.finfo(float).eps*dimension*(base_norm+abs(species_penalty)+fixed_penalty)
                    energy_envelope=64*np.finfo(float).eps*dimension*(np.linalg.norm(local_energy)+np.linalg.norm(correction))
                    if not np.isfinite(trace) or trace<=envelope or np.linalg.norm(energy_image)<=energy_envelope:
                        # Do not clip a cancellation into a positive trace. Use
                        # the original explicit global whitening/norm arithmetic.
                        whole=np.zeros((self.n,columns.shape[1]));whole[rows]=columns
                        image=whiten(whole);energy=float(np.sum((image.T@residual)**2))
                        if method=='reml':image-=left@(left.T@image)
                        trace=float(np.sum(image*image));fallbacks+=1;global_entity_cells+=whole.size
                    else:energy=float(np.sum(energy_image*energy_image))
                    traces[j]+=trace;energies[j]+=energy;local_batches+=1;local_entity_cells+=a.size
                    largest_batch=max(largest_batch,columns.shape[1])
        species_batches=0
        for start in range(0,self.factor.shape[1],self.column_batch):
            image=whiten(self.factor[:,start:start+self.column_batch])
            energies[-1]+=float(np.sum((image.T@residual)**2))
            if method=='reml':image-=left@(left.T@image)
            traces[-1]+=float(np.sum(image*image));species_batches+=1;largest_batch=max(largest_batch,image.shape[1])
        scores=.5*(traces-df*energies/quadratic)
        if not np.isfinite(scores).all():raise ArithmeticError('Unrepresentable independent score')
        return dict(method=method,negative_profiled_likelihood=float(objective),profiled_scale=scale,residual_quadratic=quadratic,
            beta=beta,conditional_beta_covariance=covariance,gradient=scores,parameter_names=self.parameter_names,
            kernel_traces_in_score=traces,residual_kernel_energies=energies,inverse_relative_residual=relative_error,
            maximum_component_condition=maximum_condition,streamed_kernel_batches=local_batches+species_batches,
            largest_dense_kernel_column_batch=largest_batch,component_local_entity_batches=local_batches,
            explicit_global_whitening_fallback_batches=fallbacks,component_local_entity_cells=local_entity_cells,
            explicit_global_entity_cells=global_entity_cells,cached_component_kernel_bytes=sum(k.nbytes for ks in self.component_kernels for k in ks if k is not None),
            coefficient_covariance_is_conditional=True,scientific_eligibility=False)
