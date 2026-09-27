#!/usr/bin/env python3
"""Check independent architecture classification across 1000 repeat/order cases."""
import random
from assess_architecture_matched_support import pair_signature
from readback_architecture_support import classify
rng=random.Random(20260927)
for _ in range(1000):
    pairs=[(tuple((f'PF{rng.randrange(3)}.1',rng.choice(['Domain','Family'])) for j in range(rng.randrange(4))),bool(rng.randrange(2))) for side in range(2)]
    assert classify(*pairs)==pair_signature(*pairs)
assert classify(((),True),((),True))==('neither_annotated',None)
print('1000 independent signature/status comparisons and empty-annotation guard passed.')
