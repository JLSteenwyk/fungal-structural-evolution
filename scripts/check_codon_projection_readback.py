#!/usr/bin/env python3
"""Cross-check independent encoded-pair intersections against direct known sets."""
import random
from readback_full_codon_projection import encoded,metrics
rng=random.Random(19)
for _ in range(300):
    universe=[(i,j) for i in range(1,8) for j in range(1,8)]
    old=set(rng.sample(universe,rng.randrange(50)));raw=set(rng.sample(universe,rng.randrange(50)));new=set(rng.sample(sorted(raw),rng.randrange(len(raw)+1)))
    def convert(s):return encoded([x for x,y in s],[y for x,y in s],8)
    r=metrics(convert(old),convert(raw),convert(new));assert r['original_preserved_raw']==len(old&raw) and r['original_preserved_retained']==len(old&new) and r['original_preserved_but_filtered']==len((old&raw)-new) and r['original_absent_from_raw']==len(old-raw) and r['local_retained_not_original']==len(new-old)
    assert r['original_pairs']==r['original_preserved_retained']+r['original_preserved_but_filtered']+r['original_absent_from_raw']
assert encoded([0,1,2],[1,0,3],8).tolist()==[19]
print('Passed 300 randomized set/encoded-array comparisons and zero-mask handling')
