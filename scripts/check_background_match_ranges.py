#!/usr/bin/env python3
"""Check graph range indexing against direct enumeration, including exact zero."""
import random
from prepare_background_match_graph import eligible_range
rng=random.Random(2718)
for _ in range(1000):
    values=sorted([0,0,1,2]+[rng.random()*10 for _ in range(40)])
    distance=rng.choice([0,1,2,rng.random()*10])
    lo,hi=eligible_range(values,distance)
    assert values[lo:hi]==[x for x in values if distance/2<=x<=distance*2]
assert eligible_range([],0)==(0,0)
for distance in [-1,float('nan'),float('inf')]:
    try:eligible_range([],distance)
    except ValueError:pass
    else:raise AssertionError('Invalid distance accepted')
print('Passed 1000 direct range enumerations, exact-zero/boundary/empty cases and invalid distances.')
