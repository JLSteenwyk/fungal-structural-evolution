#!/usr/bin/env python3
"""Known balance moments, zero/absent data handling and endpoint-invariant features."""
import numpy as np
from background_control_balance import balance,feature_vector
r=balance(np.array([1.,2.,3.]),np.array([2.,3.,4.]),np.array([0.,1.,2.,3.,4.]))
assert r['pairs']==3 and r['target_mean']==2 and r['control_mean']==3 and r['target_sd']==1 and r['control_sd']==1 and r['standardized_mean_difference']==-1 and r['mean_absolute_difference']==1 and r['selection_mean_shift']==0
assert balance(np.array([]),np.array([]),np.array([1.]))['smd_status']=='no_pairs'
assert balance(np.array([1.]),np.array([2.]),np.array([1.]))['smd_status']=='insufficient_pairs'
assert balance(np.array([1.,1.]),np.array([2.,2.]),np.array([1.,1.]))['smd_status']=='zero_pooled_variance'
r=balance(np.array([np.nan,1.,2.]),np.array([np.nan,1.,3.]),np.array([np.nan,1.,2.,3.]));assert r['pairs']==2 and r['baseline_targets']==3 and r['p95_absolute_difference']==.95
n=dict(sequence_distance=0,length_a=100,length_b=200,mean_ca_plddt_a=90,mean_ca_plddt_b=70,fraction_ca_plddt_below50_a=.1,fraction_ca_plddt_below50_b=.3)
m=dict(n)
for field in ['length','mean_ca_plddt','fraction_ca_plddt_below50']:m[field+'_a'],m[field+'_b']=n[field+'_b'],n[field+'_a']
assert np.allclose(feature_vector(n),feature_vector(m),equal_nan=True)
assert np.isnan(feature_vector(n)[1]) and feature_vector(n)[0]==0
print('Passed known moments/SMD/selection shift, absent/single/zero-variance cases, positive-log exclusion and endpoint swaps.')
