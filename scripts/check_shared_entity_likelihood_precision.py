#!/usr/bin/env python3
"""80-digit ML/REML score comparison for strongly correlated signed kernels."""
import argparse
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy import sparse
from shared_entity_likelihood import SharedEntityLikelihood
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    z=np.array([[1,0],[1,0],[-1,0],[0,1],[0,1],[0,-1]],dtype=float)
    f=np.array([[.1,.2],[.2,-.1],[.3,.3],[.4,-.2],[.5,.5],[.6,.4]])
    d=np.array([.2,.5,.7,.3,.6,.4]);x=np.column_stack([np.ones(6),np.arange(1,7)/10])
    y=np.array([.3,-.8,.1,.6,-.3,.2]);ratios=np.array([10000.,3000.])
    likelihood=SharedEntityLikelihood(['A']*3+['B']*3,{'signed':sparse.csr_matrix(z)},f,d,x,y);errors={}
    def matrix(value):return mp.matrix([[mp.mpf(str(v)) for v in row] for row in value])
    with mp.workdps(80):
        mz=matrix(z);mf=matrix(f);mx=matrix(x);my=mp.matrix([mp.mpf(str(v)) for v in y])
        kernels=[mz*mz.T,mf*mf.T];v=mp.diag([mp.mpf(str(v)) for v in d])+10000*kernels[0]+3000*kernels[1]
        iv=v**-1;info=mx.T*iv*mx;ii=info**-1;beta=ii*mx.T*iv*my
        residual=my-mx*beta;alpha=iv*residual;q=(residual.T*alpha)[0];projection=iv-iv*mx*ii*mx.T*iv
        for method in ['ml','reml']:
            df=4 if method=='reml' else 6
            objective=(mp.log(mp.det(v))+(mp.log(mp.det(info)) if method=='reml' else 0)+df*(1+mp.log(2*mp.pi*q/df)))/2
            gradients=[(sum(((projection if method=='reml' else iv)*k)[i,i] for i in range(6))-df*(alpha.T*k*alpha)[0]/q)/2 for k in kernels]
            actual=likelihood.evaluate(ratios,method)
            np.testing.assert_allclose(actual['negative_profiled_likelihood'],float(objective),rtol=3e-10,atol=3e-10)
            np.testing.assert_allclose(actual['gradient'],np.asarray(gradients,dtype=float),rtol=3e-7,atol=3e-11)
            errors[method]=dict(objective_absolute_error=abs(actual['negative_profiled_likelihood']-float(objective)),
                maximum_gradient_absolute_error=float(np.max(abs(actual['gradient']-np.asarray(gradients,dtype=float)))))
    result=dict(status='passed_80_digit_signed_shared_entity_likelihood_scores',precision_digits=80,errors=errors,
        source_hashes={str(p):sha(p) for p in [Path(__file__),Path('scripts/shared_entity_likelihood.py'),Path('scripts/shared_entity_covariance.py')]},
        scope='Strong signed numerical covariance test only; exact high-precision inverse/determinant/gradient. '
              'No production biological fit, pilot or uncertainty calibration.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
