"""Descriptive GLS marginal residual diagnostics for the matched covariance.

For a fixed covariance V, Cov(y-X beta_hat)=V-X(X'V^-1X)^-1X'.
Component estimates are fitted in the real application: these quantities are
not independent normal draws and their summaries are not calibrated tests.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from matched_mixed_covariance import MatchedCovariance


def diagnostics(background,family,factor,design,response,variances):
    x=np.asarray(design,float);y=np.asarray(response,float);f=np.asarray(factor,float)
    v=np.asarray(variances,float)
    if v.shape!=(4,):raise ValueError('Four absolute variance components required')
    c=MatchedCovariance(background,family,f,*v)
    if x.ndim!=2 or x.shape[0]!=c.n or y.shape!=(c.n,) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Incompatible or nonfinite design/response')
    n,p=x.shape
    if p==0 or n<=p or np.linalg.matrix_rank(x)!=p:raise ValueError('Full rank design and positive residual degrees of freedom required')
    inverse_x=c.solve(x);info=x.T@inverse_x;info=(info+info.T)/2
    phi=cho_solve(cho_factor(info,lower=True),np.eye(p))
    beta=phi@(x.T@c.solve(y));mean=x@beta;residual=y-mean
    diagonal=v[0]+v[1]+v[2]+v[3]*np.einsum('ij,ij->i',f,f)
    fitted_variance=np.einsum('ij,jk,ik->i',x,phi,x)
    residual_variance=diagonal-fitted_variance
    if np.any(residual_variance<=1e-12*diagonal) or not np.isfinite(residual_variance).all():
        return dict(status='residual_variance_requires_review',reason='nonpositive_or_negligible_marginal_residual_variance')
    z=residual/np.sqrt(residual_variance)
    probabilities=[0,.01,.05,.25,.5,.75,.95,.99,1]
    summaries=dict(mean=float(z.mean()),second_raw_moment=float(np.mean(z**2)),
        third_raw_moment=float(np.mean(z**3)),fourth_raw_moment=float(np.mean(z**4)),
        absolute_above_2=int(np.sum(abs(z)>2)),absolute_above_3=int(np.sum(abs(z)>3)),
        quantile_probabilities=probabilities,quantiles=np.quantile(z,probabilities).tolist(),
        maximum_fixed_effect_variance_fraction=float(np.max(fitted_variance/diagonal)))
    # Equal-count bins ordered by each covariate: descriptive heteroscedasticity/nonlinearity screens.
    bins=[]
    for j in range(p):
        if np.ptp(x[:,j])<=1e-12:continue
        order=np.argsort(x[:,j],kind='stable')
        for k,idx in enumerate(np.array_split(order,min(10,n))):
            bins.append(dict(design_column=j,bin=k,records=len(idx),covariate_min=float(x[idx,j].min()),
                covariate_max=float(x[idx,j].max()),mean_standardized_residual=float(z[idx].mean()),
                mean_squared_standardized_residual=float(np.mean(z[idx]**2))))
    return dict(status='descriptive_marginal_residual_diagnostics',records=n,coefficients=p,
        beta=beta.tolist(),summaries=summaries,covariate_bins=bins,
        scope='Marginal residuals scaled by fitted V-X Phi X prime diagonal; correlations remain. '
              'No independence, Gaussian adequacy, significance or calibrated tail probability claimed.')
