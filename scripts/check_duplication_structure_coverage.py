#!/usr/bin/env python3
"""Check coverage accounting without dropping absent models or complex events."""
import random
from inventory_duplication_structure_coverage import event_coverage

r=random.Random(20260926)
for i in range(200):
    genes=[f'F1_g{j}' for j in range(r.randrange(2,50))];cut=r.randrange(1,len(genes)+1)
    a,b=genes[:cut],genes[cut:];covered=set(r.sample(genes,r.randrange(len(genes)+1)))
    kind=r.choice(['Terminal','Non-Terminal']);row=['OG1','F1','n0','0.9',kind,', '.join(a),', '.join(b)]
    sides,found,simple,both=event_coverage(row,{g:('s','m') for g in covered})
    assert sides==[a,b] and found==[[g for g in a if g in covered],[g for g in b if g in covered]]
    assert simple==(kind=='Terminal' and len(a)==len(b)==1)
    assert both==(bool(set(a)&covered) and bool(set(b)&covered))
for row in [['x']*6,['OG1','F1','n0','1','Terminal','','F1_b'],['OG1','F1','n0','1','Terminal','F1_a','F1_a']]:
    try:event_coverage(row,{})
    except ValueError:pass
    else:raise AssertionError('Invalid event accepted')
print('Passed 200 independent coverage fixtures; rejected three malformed events.')
