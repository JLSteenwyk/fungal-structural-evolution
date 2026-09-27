#!/usr/bin/env python3
"""Exercise genetic codes, whole-codon missingness, occupancy and loss partition."""
from codon_realignment_projection import project,compare
aligned={'a':'AS-X','b':'ASDX','c':'ASDX','d':'ASDX','e':'ASDX'}
dna={t:'GCTCTG'+('GAT' if t!='a' else '')+'NNN' for t in aligned}
pos,codons,keep=project(aligned,dna,12)
assert keep==[0,1,2] and pos['a']==[1,2,0,0] and codons['a']==['GCT','CTG','???','???']
try:project(aligned,dna,1)
except AssertionError:pass
else:raise AssertionError('Wrong code accepted')
for bad in [{'a':'GCTTAA'}, {'a':'GCT'}]:
 try:project({'a':'AS'},bad,1)
 except AssertionError:pass
 else:raise AssertionError('Stop or short sequence accepted')
pos,codons,keep=project({'a':'X','b':'X'},{'a':'NNN','b':'NNN'},1);assert keep==[]
r=compare({(1,1),(2,2),(3,3)},{(1,1),(2,2),(4,4)},{(1,1),(4,4)})
assert r['original_absent_from_raw']==1 and r['original_preserved_but_filtered']==1 and r['original_preserved_retained']==1 and r['local_retained_not_original']==1
print('Passed code12/code1, gaps, ambiguous codons, exact 80-percent occupancy, empty mask, stop/length rejection and correspondence partition')
