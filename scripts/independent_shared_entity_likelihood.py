"""Independent spectral-whitening ML/REML replay, without production solvers.

Only component matrices are dense. Entity traces stream bounded column batches.
This verifies a supplied candidate, not its optimality or biological inference.
"""
import numpy as np
from scipy import sparse, linalg


class IndependentSharedEntityLikelihood:
    def __init__(self, labels, incidence, factor, diagonal, design, response, column_batch=32):
        self.labels=np.asarray(labels)
        if self.labels.ndim!=1 or not len(self.labels):raise ValueError('Nonempty labels required')
        self.n=len(self.labels)
        self.diagonal=np.asarray(diagonal,dtype=float)
        self.factor=np.asarray(factor,dtype=float)
        self.design=np.asarray(design,dtype=float)
        self.response=np.asarray(response,dtype=float)
        if not isinstance(column_batch,int) or column_batch<1:raise ValueError('Positive column batch required')
        self.column_batch=column_batch
        if self.diagonal.shape!=(self.n,) or not np.isfinite(self.diagonal).all() or np.any(self.diagonal<=0):
            raise ValueError('Finite positive residual diagonal required')
        if self.factor.ndim!=2 or self.factor.shape[0]!=self.n or not np.isfinite(self.factor).all():
            raise ValueError('Finite species factor required')
        x=self.design
        if x.ndim!=2 or x.shape[0]!=self.n or not 0<x.shape[1]<self.n or self.response.shape!=(self.n,):
            raise ValueError('Compatible active design and response required')
        if not np.isfinite(x).all() or not np.isfinite(self.response).all():raise ValueError('Finite design and response required')
        self.p=x.shape[1]
        scaled=np.asarray(x,dtype=np.longdouble)
        maximum=np.max(abs(scaled),axis=0)
        if np.any(maximum==0):raise ValueError('Inactive design columns require separate accounting')
        scaled/=maximum;length=np.sqrt(np.sum(scaled*scaled,axis=0));scaled/=length
        self.maximum=np.asarray(maximum,dtype=float);self.length=np.asarray(length,dtype=float)
        self.scaled=np.asarray(scaled,dtype=float)
        singular=linalg.svdvals(self.scaled)
        if singular[-1]<=10*max(x.shape)*np.finfo(float).eps*singular[0]:
            raise ValueError('Fixed design requires rank review')
        self.column_log_jacobian=float(2*np.sum(np.log(maximum)+np.log(length)))
        _,blocks=np.unique(self.labels,return_inverse=True)
        self.parts=[np.flatnonzero(blocks==i) for i in range(blocks.max()+1)]
        self.incidence={}
        for name,original in incidence.items():
            if name=='species':raise ValueError('Reserved species parameter name')
            z=sparse.csr_matrix(original,dtype=float,copy=True)
            if z.shape[0]!=self.n or not np.isfinite(z.data).all():raise ValueError('Invalid incidence')
            z.sum_duplicates();z.eliminate_zeros();z.sort_indices()
            # Independent column-wise test of the declared component invariant.
            columns=z.tocsc()
            for j in np.flatnonzero(np.diff(columns.indptr)):
                rows=columns.indices[columns.indptr[j]:columns.indptr[j+1]]
                if len(np.unique(blocks[rows]))!=1:raise ValueError('Entity crosses components')
            self.incidence[name]=z[:,np.unique(z.indices)].tocsr()
        self.parameter_names=[*self.incidence,'species']

    def evaluate(self, ratios, method='reml'):
        ratios=np.asarray(ratios,dtype=float)
        if ratios.shape!=(len(self.parameter_names),) or not np.isfinite(ratios).all() or np.any(ratios<0):
            raise ValueError('Explicit finite nonnegative ratios required')
        if method not in ['ml','reml']:raise ValueError('Likelihood must be ml or reml')
        roots=[];logdet=0.;maximum_component_condition=1.
        for rows in self.parts:
            base=np.diag(self.diagonal[rows])
            for ratio,z in zip(ratios[:-1],self.incidence.values()):
                if ratio:
                    local=z[rows];base+=ratio*(local@local.T).toarray()
            if not np.isfinite(base).all():raise ArithmeticError('Unrepresentable component covariance')
            values,vectors=linalg.eigh(base,driver='evd',check_finite=True)
            if values[0]<=64*np.finfo(float).eps*len(rows)*values[-1]:
                raise ArithmeticError('Component inverse square root requires precision review')
            roots.append((vectors,1/np.sqrt(values)));logdet+=float(np.log(values).sum())
            maximum_component_condition=max(maximum_component_condition,float(values[-1]/values[0]))
        def base_root(values):
            out=np.empty_like(values)
            for rows,(vectors,inverse_root) in zip(self.parts,roots):
                out[rows]=vectors@(inverse_root[:,None]*(vectors.T@values[rows]))
            return out
        if ratios[-1] and self.factor.shape[1]:
            scaled_factor=base_root(self.factor)*np.sqrt(ratios[-1])
            if not np.isfinite(scaled_factor).all():raise ArithmeticError('Unrepresentable whitened species factor')
            u,s,_=linalg.svd(scaled_factor,full_matrices=False,lapack_driver='gesvd',check_finite=True)
            inverse_species_root=1/np.hypot(1.,s)
            logdet+=float(2*np.log(np.hypot(1.,s)).sum())
        else:u=np.empty((self.n,0));inverse_species_root=np.empty(0)
        def species_root(values):
            return values+u@((inverse_species_root-1)[:,None]*(u.T@values))
        def whiten(values):return species_root(base_root(values))
        def inverse(values):return base_root(species_root(whiten(values)))
        def apply(values):
            out=self.diagonal[:,None]*values
            for ratio,z in zip(ratios[:-1],self.incidence.values()):
                if ratio:out+=ratio*(z@(z.T@values))
            if ratios[-1]:out+=ratios[-1]*(self.factor@(self.factor.T@values))
            return out
        rhs=np.column_stack([self.scaled,self.response]);inverse_rhs=inverse(rhs)
        error=float(np.max(abs(apply(inverse_rhs)-rhs),initial=0))
        relative_error=error/max(float(np.max(abs(rhs),initial=0)),np.finfo(float).tiny)
        if not np.isfinite(relative_error) or relative_error>1e-8:
            raise ArithmeticError('Independent spectral inverse failed residual check')
        wx=whiten(self.scaled);wy=whiten(self.response[:,None])[:,0]
        left,values,right=linalg.svd(wx,full_matrices=False,lapack_driver='gesvd',check_finite=True)
        if values[-1]<=10*max(wx.shape)*np.finfo(float).eps*values[0]:
            raise ArithmeticError('Whitened fixed design requires precision review')
        coefficients=right.T@((left.T@wy)/values)
        residual=wy-wx@coefficients;quadratic=float(residual@residual)
        if not np.isfinite(quadratic) or quadratic<=64*np.finfo(float).eps**2*max(float(wy@wy),np.finfo(float).tiny):
            raise ArithmeticError('Residual quadratic is zero or below numerical resolution')
        df=self.n-self.p if method=='reml' else self.n;scale=quadratic/df
        if method=='reml':logdet+=float(2*np.log(values).sum())+self.column_log_jacobian
        objective=.5*(logdet+df*(1+np.log(2*np.pi*scale)))
        beta=coefficients/self.length/self.maximum
        covariance=(right.T/(values*values))@right*scale
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            covariance=covariance/self.length[:,None]/self.maximum[:,None]/self.length[None,:]/self.maximum[None,:]
        if not np.isfinite(objective) or not np.isfinite(beta).all() or not np.isfinite(covariance).all() or np.any(np.diag(covariance)<=0):
            raise ArithmeticError('Original-unit coefficients are not representable')
        scores=[];traces=[];energies=[];batches=0;largest_batch=0
        for operator in [*self.incidence.values(),self.factor]:
            trace=energy=0.
            for start in range(0,operator.shape[1],self.column_batch):
                columns=operator[:,start:start+self.column_batch]
                columns=columns.toarray() if sparse.issparse(columns) else columns
                image=whiten(columns)
                energy+=float(np.sum((image.T@residual)**2))
                if method=='reml':image-=left@(left.T@image)
                trace+=float(np.sum(image*image));batches+=1;largest_batch=max(largest_batch,image.shape[1])
            traces.append(trace);energies.append(energy);scores.append(.5*(trace-df*energy/quadratic))
        if not np.isfinite(scores).all():raise ArithmeticError('Unrepresentable independent score')
        return dict(method=method,negative_profiled_likelihood=float(objective),profiled_scale=scale,
            residual_quadratic=quadratic,beta=beta,conditional_beta_covariance=covariance,
            gradient=np.asarray(scores),parameter_names=self.parameter_names,
            kernel_traces_in_score=np.asarray(traces),residual_kernel_energies=np.asarray(energies),
            inverse_relative_residual=relative_error,maximum_component_condition=maximum_component_condition,
            streamed_kernel_batches=batches,largest_dense_kernel_column_batch=largest_batch,
            coefficient_covariance_is_conditional=True,scientific_eligibility=False)


def audit_candidate(independent, candidate, objective_atol=1e-7, gradient_atol=1e-6):
    """Replay a candidate; no optimization, curvature or uncertainty acceptance."""
    if not isinstance(independent,IndependentSharedEntityLikelihood):raise TypeError('Independent likelihood required')
    if not np.isfinite(objective_atol) or objective_atol<=0 or not np.isfinite(gradient_atol) or gradient_atol<=0:
        raise ValueError('Positive finite comparison tolerances required')
    if candidate['parameter_names']!=independent.parameter_names:raise ValueError('Candidate parameter identity mismatch')
    if candidate.get('scientific_eligibility') is not False or candidate.get('coefficient_covariance_is_conditional') is not True:
        raise ValueError('Numerical candidate inference/conditional covariance flags invalid')
    ratios=np.asarray(candidate['variance_ratios'],dtype=float)
    replay=independent.evaluate(ratios,candidate['method'])
    objective=float(candidate['negative_profiled_likelihood'])
    if not np.isfinite(objective):raise ValueError('Nonfinite candidate likelihood')
    objective_error=abs(replay['negative_profiled_likelihood']-objective)
    covariance=np.asarray(candidate['conditional_beta_covariance'])
    if covariance.shape!=(independent.p,independent.p) or not np.isfinite(covariance).all():raise ValueError('Invalid candidate coefficient covariance')
    beta=np.asarray(candidate['beta'])
    if beta.shape!=(independent.p,) or not np.isfinite(beta).all():raise ValueError('Invalid candidate coefficients')
    beta_error=float(np.max(abs(replay['beta']-beta),initial=0))
    # Compare in balanced design units. Absolute original-unit tolerances can
    # silently accept gross relative corruption of tiny-unit coefficients.
    def balanced_beta(value):return value*independent.maximum*independent.length
    def balanced_covariance(value):
        return (value*independent.maximum[:,None]*independent.length[:,None]
            *independent.maximum[None,:]*independent.length[None,:])/replay['profiled_scale']
    balanced_error=float(np.linalg.norm(balanced_beta(replay['beta'])-balanced_beta(beta)))
    coefficient_tolerance=1e-7*max(float(np.linalg.norm(balanced_beta(replay['beta']))),
        float(np.sqrt(replay['profiled_scale'])),np.finfo(float).tiny)
    coefficients_match=(np.isfinite(balanced_error) and balanced_error<=coefficient_tolerance and
        np.allclose(balanced_covariance(replay['conditional_beta_covariance']),balanced_covariance(covariance),rtol=1e-6,atol=1e-10))
    scale=float(candidate['profiled_scale']);components=np.asarray(candidate['variance_components'],dtype=float)
    if not np.isfinite(scale) or scale<=0 or components.shape!=ratios.shape or not np.isfinite(components).all() or np.any(components<0):
        raise ValueError('Invalid candidate scale/components')
    scale_match=bool(np.isclose(scale,replay['profiled_scale'],rtol=1e-7,atol=0))
    components_match=bool(np.allclose(components,replay['profiled_scale']*ratios,rtol=1e-7,atol=0))
    norms=np.asarray(candidate['kernel_normalization'],dtype=float)
    expected_norms=np.asarray([float(z.multiply(z).sum())/independent.n for z in independent.incidence.values()]
        +[float(np.sum(independent.factor*independent.factor))/independent.n])/float(np.mean(independent.diagonal))
    if norms.shape!=ratios.shape or not np.isfinite(norms).all() or np.any(norms<=0) or not np.allclose(norms,expected_norms,rtol=1e-12,atol=0):
        raise ValueError('Candidate kernel normalization mismatch')
    coordinate_gradient=replay['gradient']*(1+ratios*norms)/norms
    projected=coordinate_gradient.copy();projected[ratios==0]=np.minimum(projected[ratios==0],0)
    kkt=float(np.max(abs(projected),initial=0))
    cap=float(candidate['maximum_scaled_variance'])
    if not np.isfinite(cap) or cap<=0:raise ValueError('Candidate variance cap invalid')
    upper_hit=bool(np.any(ratios*norms>=cap*(1-64*np.finfo(float).eps)))
    passed=(np.isfinite(objective_error) and objective_error<=objective_atol and coefficients_match and scale_match and components_match and
        kkt<=gradient_atol and not upper_hit and candidate['status']=='optimized_shared_entity_candidate_pending_independent_audit')
    return dict(status='independent_shared_entity_replay_passed_pending_curvature_and_calibration' if passed
        else 'independent_shared_entity_replay_requires_review',objective_absolute_error=float(objective_error),
        maximum_coefficient_absolute_error=beta_error,coefficient_and_conditional_covariance_comparison_passed=bool(coefficients_match),
        balanced_coefficient_error_norm=balanced_error,balanced_coefficient_tolerance=coefficient_tolerance,
        profiled_scale_comparison_passed=scale_match,variance_components_comparison_passed=components_match,
        maximum_independent_projected_gradient=kkt,upper_variance_cap_hit=upper_hit,
        inverse_relative_residual=replay['inverse_relative_residual'],objective_absolute_tolerance=objective_atol,
        gradient_absolute_tolerance=gradient_atol,scientific_eligibility=False,
        scope='Independent spectral likelihood/conditional coefficients and variance boundary KKT replay only. '
            'No independent global optimum, curvature, uncertainty calibration or accepted biological effects.')


def audit_local_curvature(independent,candidate,coordinate_step=1e-4,gradient_atol=1e-6):
    """Two-resolution score Hessian, with a sufficient constrained-minimum check.

    Weak lower boundaries stay in the tested space. Strictly positive outward
    scores permit removing their coordinates from the critical cone. Positive
    curvature on the whole remaining space is sufficient, not necessary.
    """
    if not np.isfinite(coordinate_step) or coordinate_step<=0:raise ValueError('Positive finite curvature step required')
    replay=audit_candidate(independent,candidate,gradient_atol=gradient_atol)
    if replay['status']!='independent_shared_entity_replay_passed_pending_curvature_and_calibration':
        return dict(status='independent_shared_entity_curvature_requires_review',candidate_replay=replay,
            curvature_evaluations=0,scientific_eligibility=False,reason='Candidate failed independent replay')
    ratios=np.asarray(candidate['variance_ratios']);norms=np.asarray(candidate['kernel_normalization'])
    coordinates=np.log1p(ratios*norms);upper=np.log1p(candidate['maximum_scaled_variance']);evaluations=0
    def score(point):
        nonlocal evaluations
        if np.any(point<0) or np.any(point>upper):raise ArithmeticError('Curvature point outside explicit variance bounds')
        result=independent.evaluate(np.expm1(point)/norms,candidate['method']);evaluations+=1
        return result['gradient']*np.exp(point)/norms
    initial=score(coordinates);matrices=[];stencils=[];minimum_step=float('inf');score_max=float(np.max(abs(initial),initial=0))
    for resolution in [1.,.5]:
        hessian=np.empty((len(norms),len(norms)));used=[]
        for j,value in enumerate(coordinates):
            step=coordinate_step*resolution*max(1.,value)
            # Two one-sided steps must fit even in an explicit narrow box.
            step=min(step,upper/4)
            if step<=128*np.finfo(float).eps*max(1.,upper):raise ArithmeticError('Curvature finite difference below resolution')
            minimum_step=min(minimum_step,step)
            def shifted(distance):
                point=coordinates.copy();point[j]+=distance*step;result=score(point)
                nonlocal score_max
                score_max=max(score_max,float(np.max(abs(result),initial=0)))
                return result
            if value>=step and upper-value>=step:
                hessian[:,j]=(shifted(1)-shifted(-1))/(2*step);used.append('central')
            elif upper-value>=2*step:
                hessian[:,j]=(-3*initial+4*shifted(1)-shifted(2))/(2*step);used.append('forward')
            elif value>=2*step:
                hessian[:,j]=(3*initial-4*shifted(-1)+shifted(-2))/(2*step);used.append('backward')
            else:raise ArithmeticError('No resolvable bounded curvature stencil')
        if not np.isfinite(hessian).all():raise ArithmeticError('Nonfinite independent curvature')
        matrices.append(hessian);stencils.append(used)
    antisymmetry=float(np.linalg.norm(matrices[-1]-matrices[-1].T,ord=2))
    discrepancy=float(np.linalg.norm(matrices[-1]-matrices[0],ord=2))
    numerical_envelope=4*discrepancy+128*np.finfo(float).eps*max(1.,score_max)*len(norms)/minimum_step
    strong_lower=(coordinates==0)&(initial>gradient_atol)
    tested=np.flatnonzero(~strong_lower);matrix=(matrices[-1]+matrices[-1].T)/2
    eigenvalues=linalg.eigvalsh(matrix[np.ix_(tested,tested)]) if len(tested) else np.empty(0)
    qualified=(antisymmetry<=numerical_envelope and (not len(eigenvalues) or eigenvalues[0]>numerical_envelope))
    return dict(status='independent_shared_entity_local_minimum_check_passed_pending_global_and_inferential_audits' if qualified
        else 'independent_shared_entity_curvature_requires_review',candidate_replay=replay,
        curvature_evaluations=evaluations,coordinate_step=coordinate_step,finite_difference_stencils=stencils,
        score_hessian=matrix.tolist(),strong_lower_boundary_indices=np.flatnonzero(strong_lower).tolist(),
        tested_critical_space_indices=tested.tolist(),tested_curvature_eigenvalues=eigenvalues.tolist(),
        hessian_antisymmetry_norm=antisymmetry,two_resolution_hessian_difference_norm=discrepancy,
        numerical_curvature_review_envelope=float(numerical_envelope),scientific_eligibility=False,
        scope='Sufficient local constrained minimum check from independent spectral scores at two '
            'finite-difference resolutions. Positive curvature required on a superset of the critical '
            'cone; flat/indefinite/asymmetric/unresolved cases remain review. Numerical envelope is '
            'an empirical stability screen, not a rigorous floating-point error bound. No global '
            'optimum, variance uncertainty, likelihood test calibration or biological acceptance.')
