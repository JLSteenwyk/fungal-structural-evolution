#!/usr/bin/env python3
"""Profile ordinary Gaussian ML for comparing matched fixed-effect designs."""
import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from matched_mixed_covariance import MatchedCovariance, profiled_reml as _direct_reml


def profiled_ml(covariance, design, response):
    """Use the direct residual path; beta is identical at fixed covariance ratios."""
    result = _direct_reml(covariance, design, response)
    q = result['residual_quadratic']
    old_scale = result['profiled_scale']
    scale = q / covariance.n
    result.pop('negative_profiled_reml')
    result['negative_profiled_ml'] = float(.5 * (covariance.logdet + covariance.n * (1 + np.log(2*np.pi*scale))))
    result['profiled_scale'] = float(scale)
    result['conditional_beta_covariance'] *= scale / old_scale
    result['variance_profile_denominator'] = covariance.n
    return result


class CachedMatchedML:
    def __init__(self, background, family, factor, design, response):
        # Preserve copies for a numerically stable fallback near exact fits.
        self.background=np.asarray(background).copy();self.family=np.asarray(family).copy()
        self.factor=np.asarray(factor,dtype=float).copy()
        self.x=np.asarray(design,dtype=float).copy();self.y=np.asarray(response,dtype=float).copy()
        validation=MatchedCovariance(self.background,self.family,self.factor)
        self.n=validation.n
        if self.x.ndim!=2 or self.x.shape[0]!=self.n or self.y.shape!=(self.n,):
            raise ValueError('Incompatible fixed design or response')
        self.p=self.x.shape[1];self.r=self.factor.shape[1]
        if not np.isfinite(self.x).all() or not np.isfinite(self.y).all() or self.p==0 or self.n<=self.p:
            raise ValueError('Finite full-rank design and positive residual degrees of freedom required')
        if np.linalg.matrix_rank(self.x)!=self.p:raise ValueError('Rank-deficient fixed design')
        bg,fam=validation.bg,validation.fam
        self.counts=np.bincount(bg)
        group_family=np.empty(validation.ng,dtype=int)
        group_family[bg]=fam
        self.family_incidence=sparse.csr_matrix((np.ones(validation.ng),(group_family,np.arange(validation.ng))),shape=(validation.nf,validation.ng))
        matrix=np.column_stack([self.factor,self.x,self.y])
        self.gram=matrix.T@matrix
        self.sums=np.zeros((validation.ng,matrix.shape[1]))
        np.add.at(self.sums,bg,matrix)
        self.distinct_counts=np.unique(self.counts)
        self.group_grams=np.array([self.sums[self.counts==c].T@self.sums[self.counts==c] for c in self.distinct_counts])
        self.group_count_multiplicities=np.array([np.sum(self.counts==c) for c in self.distinct_counts])

    def evaluate(self, background_variance, family_variance, species_variance):
        b,f,p=map(float,[background_variance,family_variance,species_variance])
        if not np.isfinite([b,f,p]).all() or min(b,f,p)<0:raise ValueError('Invalid variance ratios')
        offset=0 if p>0 and self.r else self.r
        weights=1./(1.+b*self.counts)
        coefficients=b/(1.+b*self.distinct_counts)
        gram=self.gram[offset:,offset:].copy()
        if b:
            gram-=np.tensordot(coefficients,self.group_grams[:,offset:,offset:],axes=(0,0))
        logdet=float(np.sum(self.group_count_multiplicities*np.log1p(b*self.distinct_counts)))
        if f:
            family_sums=self.family_incidence@(weights[:,None]*self.sums[:,offset:])
            information=self.family_incidence@(weights*self.counts)
            family_coefficient=f/(1.+f*information)
            gram-=family_sums.T@(family_coefficient[:,None]*family_sums)
            logdet+=float(np.log1p(f*information).sum())
        if offset==0 and p>0 and self.r:
            core=np.eye(self.r)+p*gram[:self.r,:self.r]
            core=(core+core.T)/2
            chol=cho_factor(core,lower=True)
            cross=gram[:self.r,self.r:]
            gram=gram[self.r:,self.r:]-p*cross.T@cho_solve(chol,cross)
            logdet+=float(2*np.log(np.diag(chol[0])).sum())
        gram=(gram+gram.T)/2
        if not np.isfinite(gram).all() or not np.isfinite(logdet):
            raise FloatingPointError('Nonfinite covariance cross-products or determinant')
        information=gram[:self.p,:self.p]
        chol=cho_factor(information,lower=True)
        xy=gram[:self.p,-1]
        beta=cho_solve(chol,xy)
        yy=float(gram[-1,-1]);explained=float(beta@xy)
        quadratic=yy-explained
        if not np.isfinite(quadratic):raise FloatingPointError('Nonfinite residual quadratic form')
        # Cross-product subtraction loses relative precision near a perfect fit.
        if quadratic<=1e-9*max(abs(yy)+abs(explained),np.finfo(float).tiny):
            result=profiled_ml(MatchedCovariance(self.background,self.family,self.factor,1.,b,f,p),self.x,self.y)
            result['evaluation_method']='direct_residual_fallback'
            return result
        scale=quadratic/self.n
        objective=.5*(logdet+self.n*(1+np.log(2*np.pi*scale)))
        return dict(negative_profiled_ml=float(objective),beta=beta,profiled_scale=float(scale),
                    residual_quadratic=float(quadratic),conditional_beta_covariance=scale*cho_solve(chol,np.eye(self.p)),
                    residual_degrees_of_freedom=self.n-self.p,variance_profile_denominator=self.n,evaluation_method='cached_cross_products')
