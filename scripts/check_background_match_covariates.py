#!/usr/bin/env python3
"""Check endpoint swaps, joint-orientation requirements, overlaps and invalid values."""
from background_match_covariates import differences

def node(lengths,confidence):
    n={}
    for s,l,c in zip('ab',lengths,confidence):
        n.update({k+'_'+s:v for k,v in dict(length=l,mean_ca_plddt=c,fraction_ca_plddt_below50=0,gene='g'+s,model_id='m'+s,version=1,sequence_sha256='s'+s).items()})
    return n
x=node([100,200],[90,50]);r=differences(x,x);assert all(r['within_'+b]==1 for b in ['tight','moderate','wide']) and r['shared_genes']==r['shared_models']==r['shared_sequences']==2
# One orientation matches lengths and the other confidence: neither is jointly acceptable.
y=node([100,200],[50,90]);r=differences(x,y);assert r['length_ratio_0']==1 and r['plddt_difference_1']==0 and r['within_wide']==0
z={k[:-1]+('b' if k[-1]=='a' else 'a'):v for k,v in y.items()};s=differences(x,z)
for field in ['length_ratio','plddt_difference','lowconf_difference']:assert r[field+'_0']==s[field+'_1'] and r[field+'_1']==s[field+'_0']
for band in ['tight','moderate','wide']:assert r['within_'+band]==s['within_'+band]
for key,value in [('length_a',0),('mean_ca_plddt_a',101),('fraction_ca_plddt_below50_b',-1),('length_b',float('nan'))]:
    bad=dict(x);bad[key]=value
    try:differences(x,bad)
    except ValueError:pass
    else:raise AssertionError('Invalid covariate accepted')
print('Passed identity, endpoint permutation, joint-orientation counterexample, overlap and invalid-value checks.')
