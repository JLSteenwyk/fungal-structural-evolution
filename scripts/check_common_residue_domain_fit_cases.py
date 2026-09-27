#!/usr/bin/env python3
"""Check common-core distances, duplicate-label orientation and rigid-motion invariance."""
import numpy as np
from fit_domain_triad_common_residues import fit_triplet
rng=np.random.default_rng(4131);r=rng.normal(size=(50,3));a=r+rng.normal(size=r.shape)*.1;b=r+rng.normal(size=r.shape)*.4
letters=['A'*50,'A'*25+'C'*25,'C'*50];confidence=[np.full(50,90.),np.full(50,80.),np.full(50,60.)]
x=fit_triplet([a,b,r],letters,confidence)
assert x['fit_status']=='computed_unique_at_numeric_tolerance'
assert abs(x['rmsd_ar_minus_br'])<=x['rmsd_ab']+1e-12
assert x['sequence_identity_ab']==.5 and x['sequence_identity_ar']==0 and x['sequence_identity_br']==.5
assert x['joint_plddt70_fraction']==0
s=fit_triplet([b,a,r],[letters[1],letters[0],letters[2]],[confidence[1],confidence[0],confidence[2]])
assert np.isclose(x['rmsd_ar_minus_br'],-s['rmsd_ar_minus_br'],atol=1e-12)
assert np.isclose(x['rmsd_ab'],s['rmsd_ab'],atol=1e-12)
q,_=np.linalg.qr(rng.normal(size=(3,3)));q[:,0]*=np.linalg.det(q)
t=fit_triplet([a@q+[5,1,9],b+[-9,2,3],r],letters,confidence)
for k in ['rmsd_ab','rmsd_ar','rmsd_br','rmsd_ar_minus_br']:assert np.isclose(x[k],t[k],atol=1e-12)
line=np.column_stack([np.arange(50),np.zeros(50),np.zeros(50)])
assert fit_triplet([line]*3,letters,confidence)['fit_status']=='computed_degenerate_geometry'
try:fit_triplet([r[:2]]*3,['AA']*3,[np.ones(2)]*3)
except ValueError:pass
else:raise AssertionError('Two-residue core accepted')
print('Same-core metric bounds, label reversal, rigid motion, sequence identity, confidence and short/degenerate-core checks passed.')
