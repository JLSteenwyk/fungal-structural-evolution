#!/usr/bin/env python3
"""Analytic profiled-REML gradient for diagnostics; does not change live fits."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def evaluate_gradient(cache, ratios):
    b,f,p=np.asarray(ratios,dtype=float)
    if not np.isfinite([b,f,p]).all() or min(b,f,p)<0:raise ValueError('Invalid ratios')
    n=cache.counts; sizes=cache.distinct_counts
    w=1/(1+b*n); dw=-n*w*w
    t=cache.family_incidence@(w[:,None]*cache.sums)
    tb=cache.family_incidence@(dw[:,None]*cache.sums)
    d=cache.family_incidence@(n*w); db=cache.family_incidence@(n*dw)
    denom=1+f*d; alpha=f/denom
    alphab=-f*f*db/(denom*denom);alphaf=1/(denom*denom)
    gram=cache.gram-np.tensordot(b/(1+b*sizes),cache.group_grams,axes=(0,0))-t.T@(alpha[:,None]*t)
    gb=-np.tensordot(1/(1+b*sizes)**2,cache.group_grams,axes=(0,0))-tb.T@(alpha[:,None]*t)-t.T@(alpha[:,None]*tb)-t.T@(alphab[:,None]*t)
    gf=-t.T@(alphaf[:,None]*t)
    logdet=float(np.sum(cache.group_count_multiplicities*np.log1p(b*sizes))+np.log1p(f*d).sum())
    ldb=float(np.sum(cache.group_count_multiplicities*sizes/(1+b*sizes))+np.sum(f*db/denom))
    ldf=float(np.sum(d/denom))
    r=cache.r
    k=gram[:r,:r];cross=gram[:r,r:]
    if r:
        core=np.eye(r)+p*k;core=(core+core.T)/2
        ch=cho_factor(core,lower=True);inverse=cho_solve(ch,np.eye(r));u=cho_solve(ch,cross)
        logdet+=float(2*np.log(np.diag(ch[0])).sum())
    else:
        inverse=np.zeros((0,0));u=np.zeros_like(cross)
    reduced=gram[r:,r:]-p*cross.T@u
    derivatives=[]
    logderivatives=[]
    for derivative,ld in [(gb,ldb),(gf,ldf)]:
        v=derivative[:r,r:]
        derivatives.append(derivative[r:,r:]-p*(v.T@u+u.T@v)+p*p*u.T@derivative[:r,:r]@u)
        logderivatives.append(ld+p*float(np.sum(inverse*derivative[:r,:r].T)))
    derivatives.append(-u.T@u)
    logderivatives.append(float(np.sum(inverse*k.T)))
    reduced=(reduced+reduced.T)/2
    information=reduced[:cache.p,:cache.p]
    ch=cho_factor(information,lower=True);inv_info=cho_solve(ch,np.eye(cache.p))
    beta=cho_solve(ch,reduced[:cache.p,-1])
    explained=float(beta@reduced[:cache.p,-1]);yy=float(reduced[-1,-1]);q=yy-explained
    if not np.isfinite(q) or q<=1e-9*max(abs(yy)+abs(explained),np.finfo(float).tiny):
        raise ArithmeticError('Cross-product cancellation: analytic diagnostic requires separate residual treatment')
    residual_contrast=np.r_[-beta,1.];df=cache.n-cache.p
    gradient=np.array([.5*(ld+np.sum(inv_info*dg[:cache.p,:cache.p].T)+df*float(residual_contrast@dg@residual_contrast)/q)
                       for dg,ld in zip(derivatives,logderivatives)])
    objective=.5*(logdet+2*np.log(np.diag(ch[0])).sum()+df*(1+np.log(2*np.pi*q/df)))
    if not np.isfinite(gradient).all() or not np.isfinite(objective):raise FloatingPointError('Nonfinite analytic gradient')
    return dict(negative_profiled_reml=float(objective),ratio_gradient=gradient,
                log1p_ratio_gradient=gradient*(1+np.asarray(ratios)))
