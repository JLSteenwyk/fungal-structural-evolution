#!/usr/bin/env python3
"""Binary gap likelihood with direct probability of an observed-one event."""
import math
import numpy as np
from scipy.special import gammainc, logsumexp
from scipy.stats import gamma
from fit_conditional_indel_models import Likelihood


class StableLikelihood(Likelihood):
    def __init__(self,tree,sequences,correction):
        super().__init__(tree,sequences,correction)
        masks=(self.patterns=='?') if correction=='observed_mask' else np.zeros(self.patterns.shape,dtype=bool)
        self.masks,self.mask_inverse=np.unique(masks,axis=0,return_inverse=True)
        self.tip_column={self.index[n]:i for i,n in enumerate(tree.get_terminals())}

    def evaluate(self,parameters):
        gain,loss,alpha=np.exp(parameters)
        pi=np.array([loss,gain])/(gain+loss)
        rates=4*np.diff(gammainc(alpha+1,alpha*gamma.ppf([0,.25,.5,.75,1],a=alpha,scale=1/alpha)))
        observed_logs,event_logs=[],[]
        for rate in rates:
            values,scales,absent,present={},{},{},{}
            for i,children in enumerate(self.children):
                part=self.evidence[i].copy() if i in self.evidence else np.ones((2,self.width))
                scale=np.zeros(self.width)
                no=np.ones((2,len(self.masks)))
                yes=np.zeros_like(no)
                if i in self.tip_column:
                    mask=self.masks[:,self.tip_column[i]]
                    no[1,~mask]=0.
                    yes[1,~mask]=1.
                for child in children:
                    changed=-np.expm1(-(gain+loss)*rate*self.lengths[child])
                    transition=np.array([[1-pi[1]*changed,pi[1]*changed],
                                         [pi[0]*changed,1-pi[0]*changed]])
                    part*=transition@values.pop(child)
                    scale+=scales.pop(child)
                    norm=part.max(axis=0)
                    if np.any(norm<=0):
                        return -np.inf
                    part/=norm
                    scale+=np.log(norm)
                    child_no=transition@absent.pop(child)
                    child_yes=transition@present.pop(child)
                    # Union of disjoint events, without subtracting near-unity values.
                    yes=yes*(child_no+child_yes)+no*child_yes
                    no=no*child_no
                values[i],scales[i]=part,scale
                absent[i],present[i]=no,yes
            root=len(self.nodes)-1
            observed_logs.append(np.log(pi@values[root])+scales[root])
            event=pi@present[root]
            assert np.all(event>0) and np.all(event<=1+1e-12)
            event_logs.append(np.log(event))
        log_numerator=logsumexp(observed_logs,axis=0)-math.log(4)
        log_denominator=logsumexp(event_logs,axis=0)-math.log(4)
        return float(self.counts@(log_numerator[self.observed]-log_denominator[self.mask_inverse]))
