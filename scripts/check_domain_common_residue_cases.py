#!/usr/bin/env python3
"""Test reference intersections, transitive correspondence and nonbijective rejection."""
from prepare_domain_triad_common_residues import common_residues

ab={10:20,11:21,12:22};ar={100:10,101:11,102:12};br={100:20,101:21,103:23}
common,cycle=common_residues(ab,ar,br)
assert common==cycle==[[10,20,100],[11,21,101]]
common,cycle=common_residues({10:21,11:20},ar,br)
assert len(common)==2 and cycle==[]
assert common_residues({},ar,{200:20})==([],[])
assert common_residues({},ar,br)[1]==[]
# Original positions need not be contiguous and must not be replaced by indices.
assert common_residues({901:702},{3000:901},{3000:702})==([[901,702,3000]],[[901,702,3000]])
for maps in [({1:2,3:2},{},{ }),({},{1:2,3:2},{}),({},{},{1:2,3:2})]:
    try:common_residues(*maps)
    except ValueError:pass
    else:raise AssertionError('Nonbijective mapping accepted')
# Swapping duplicate labels swaps the first two columns, not reference positions.
c,z=common_residues({v:k for k,v in ab.items()},br,ar)
assert c==z==[[20,10,100],[21,11,101]]
print('Common-reference, inconsistent-cycle, disjoint/missing maps, original-position and duplicate-label tests passed.')
