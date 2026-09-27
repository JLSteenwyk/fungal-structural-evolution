#!/usr/bin/env python3
"""Exact conditioning of a Bernoulli working distribution on compatible gap runs.

This changes supplied marginals. It is not a fitted joint evolutionary indel model.
"""
from bisect import bisect_left
import numpy as np


class CompatibleGaps:
    def __init__(self,intervals,probabilities):
        p=np.asarray(probabilities,dtype=float)
        assert len(intervals)==len(p) and np.isfinite(p).all() and np.all((p>=0)&(p<=1))
        assert all(a<=b for a,b in intervals) and len(set(intervals))==len(intervals)
        self.order=np.array(sorted(range(len(p)),key=lambda i:(intervals[i][1],intervals[i][0])),dtype=int)
        self.intervals=[intervals[i] for i in self.order];p=p[self.order];n=len(p)
        with np.errstate(divide='ignore'):
            lp=np.log(p);lq=np.log1p(-p)
        zeros=np.r_[0,np.cumsum(np.isneginf(lq))]
        finite=np.r_[0,np.cumsum(np.where(np.isfinite(lq),lq,0))]
        self.forward=np.full(n+1,-np.inf);self.forward[0]=0
        self.maximum=np.full(n+1,-np.inf);self.maximum[0]=0
        self.previous=np.zeros(n,dtype=int);self.take=np.full(n,-np.inf);self.skip=lq
        ends=[b for a,b in self.intervals]
        for j,(a,b) in enumerate(self.intervals):
            previous=bisect_left(ends,a-1,0,j)
            self.previous[j]=previous
            excluded=-np.inf if zeros[j]>zeros[previous] else finite[j]-finite[previous]
            self.take[j]=lp[j]+excluded
            self.forward[j+1]=np.logaddexp(lq[j]+self.forward[j],self.take[j]+self.forward[previous])
            self.maximum[j+1]=max(lq[j]+self.maximum[j],self.take[j]+self.maximum[previous])
        self.log_compatibility_probability=float(self.forward[-1])
        self.feasible=bool(np.isfinite(self.forward[-1]))

    def marginals(self):
        if not self.feasible:raise ValueError('zero probability of any compatible configuration')
        n=len(self.order);mass=np.zeros(n+1);mass[n]=1.;selected=np.zeros(n)
        for k in range(n,0,-1):
            if mass[k]==0:continue
            j=k-1;previous=self.previous[j]
            take=np.exp(self.take[j]+self.forward[previous]-self.forward[k])
            skip=np.exp(self.skip[j]+self.forward[j]-self.forward[k])
            assert abs(take+skip-1)<1e-10
            selected[j]=mass[k]*take
            mass[previous]+=mass[k]*take;mass[j]+=mass[k]*skip
        assert abs(mass[0]-1)<1e-10
        result=np.empty(n);result[self.order]=selected
        return result

    def configuration(self,rng=None):
        """MAP if rng=None, otherwise an exact draw from the conditioned surrogate."""
        if not self.feasible:raise ValueError('zero probability of any compatible configuration')
        chosen=[];k=len(self.order)
        while k:
            j=k-1;previous=self.previous[j]
            if rng is None:
                take=self.take[j]+self.maximum[previous]>self.skip[j]+self.maximum[j]
            else:
                take=rng.random()<np.exp(self.take[j]+self.forward[previous]-self.forward[k])
            if take:chosen.append(int(self.order[j]));k=previous
            else:k=j
        return sorted(chosen)
