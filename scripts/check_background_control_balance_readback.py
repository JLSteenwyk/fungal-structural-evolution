#!/usr/bin/env python3
"""Cross-check independent balance reconstruction on random and nonestimable cases."""
import math
import numpy as np
from background_control_balance import balance
from readback_background_control_balance import reconstruct
rng=np.random.default_rng(92728)
cases=[(np.array([]),np.array([]),np.array([])),(np.array([1.]),np.array([2.]),np.array([1.])),(np.array([1.,1.]),np.array([2.,2.]),np.array([1.,1.]))]
for i in range(250):
    x=rng.normal(size=i+2);y=rng.normal(size=i+2);base=rng.normal(size=i+5)
    if i%3==0:x[0]=y[0]=base[0]=np.nan
    cases.append((x,y,base))
for x,y,b in cases:
    a=balance(x,y,b);r=reconstruct(x,y,b);assert set(a)==set(r)
    for k,v in a.items():
        if isinstance(v,float):assert math.isclose(v,r[k],rel_tol=1e-9,abs_tol=1e-10),(k,v,r[k])
        else:assert v==r[k],(k,v,r[k])
print('Passed253 full balance reconstructions including absent/single/constant/NaN cases.')
