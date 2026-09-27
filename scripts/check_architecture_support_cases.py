#!/usr/bin/env python3
"""Check conservative whole-architecture support labels without dropping unknowns."""
from assess_architecture_matched_support import pair_signature
A=(('PF1.1','Domain'),);B=(('PF2.1','Domain'),)
assert pair_signature((A,True),(A,True))==('shared_conservative_architecture',A)
for left,right,status in [(((),False),((),False),'neither_annotated'),((A,True),((),False),'one_unannotated'),((A,True),(B,True),'different_ordered_annotations'),((A,True),(A,False),'shared_but_nonconservative'),((A+A,True),(A,True),'different_ordered_annotations'),((A+B,True),(B+A,True),'different_ordered_annotations')]:
    assert pair_signature(left,right)==(status,None)
print('Shared, missing, nonconservative, content, repeat and order cases passed.')
