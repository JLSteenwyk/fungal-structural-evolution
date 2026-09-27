#!/usr/bin/env python3
"""Compare independent array reconstruction with scalar records under random permutations."""
import numpy as np
from background_match_covariates import differences
from readback_background_match_covariates import matrix_differences
rng=np.random.default_rng(927)
t=rng.uniform(size=(1000,2,3));b=rng.uniform(size=(1000,2,3))
for a in [t,b]:a[:,:,0]=a[:,:,0]*1000+1;a[:,:,1]*=100
result=matrix_differences(t,b)
for i in range(1000):
    nodes=[]
    for a in [t[i],b[i]]:
        n={}
        for side,v in zip('ab',a):
            for field,x in zip(['length','mean_ca_plddt','fraction_ca_plddt_below50'],v):n[field+'_'+side]=float(x)
            n.update({'gene_'+side:side,'model_id_'+side:side,'version_'+side:1,'sequence_sha256_'+side:side})
        nodes.append(n)
    expected=differences(*nodes)
    for field,values in result.items():assert np.isclose(values[i],expected[field],rtol=1e-13,atol=1e-13)
swap=matrix_differences(t,b[:,::-1,:])
for band in ['tight','moderate','wide']:assert np.array_equal(result['within_'+band],swap['within_'+band])
print('Passed 1000 scalar/independent-array comparisons and endpoint-swap invariance.')
