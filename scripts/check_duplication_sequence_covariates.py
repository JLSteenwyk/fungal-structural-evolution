#!/usr/bin/env python3
"""Known distances, scale/cancellation and invalid-source checks."""
from prepare_duplication_sequence_covariates import derive
import math


def row(a,b,r):
    return dict(distance_a_to_reference=a+r, distance_b_to_reference=b+r,
                duplicate_pair_sequence_distance=a+b, reference_distance_from_duplicate_node=r)

x=derive(row(1,3,7))
assert x['sequence_tip_a']==1 and x['sequence_tip_b']==3
assert x['sequence_normalized_contrast']==-.5 and x['sequence_direction']=='b_longer'
y=derive(row(3,1,7));assert y['sequence_normalized_contrast']==.5
z=derive(row(0,0,100));assert z['sequence_normalized_contrast']==''
assert derive(row(1e-14,2e-14,100))['sequence_direction']=='unresolved'
assert derive(row(2,2,0))['sequence_direction']=='unresolved'
for scale in [.001,1,1000]:
    assert math.isclose(derive(row(scale,3*scale,7*scale))['sequence_normalized_contrast'],-.5)
for key,value in [('distance_a_to_reference',-1),('distance_b_to_reference',float('nan')),('duplicate_pair_sequence_distance',100),('reference_distance_from_duplicate_node',99)]:
    bad=row(1,3,7);bad[key]=value
    try:derive(bad)
    except ValueError:pass
    else:raise AssertionError('Bad source accepted')
print('Known contrasts, reversal, scale, near-zero, equality and invalid paths passed.')
