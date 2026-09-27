#!/usr/bin/env python3
"""Check endpoint normalization under unequal domain lengths and reversed order."""
from prepare_matched_domain_contrasts import measurement
n=dict(aligned_length='40',rmsd_recomputed='1.2',sequence_identity_exact='.7',joint_plddt70_fraction='.8',tm_left_native='.3',tm_right_native='.6',coverage_left='.4',coverage_right='.8')
c={'order0_original_coverage_a':'.2','order0_original_coverage_b':'.4','order1_original_coverage_a':'.3','order1_original_coverage_b':'.6'}
a=measurement(n,c,0);b=measurement(n,c,1)
assert (a['tm_a'],a['tm_b'],a['retained_coverage_a'],a['retained_coverage_b'])==(.3,.6,.4,.8)
assert (b['tm_a'],b['tm_b'],b['retained_coverage_a'],b['retained_coverage_b'])==(.6,.3,.8,.4)
assert (a['original_coverage_a'],a['original_coverage_b'])==(.2,.4)
assert (b['original_coverage_a'],b['original_coverage_b'])==(.3,.6)
assert a['rmsd_recomputed']==b['rmsd_recomputed']==1.2 and a['aligned_length']==b['aligned_length']==40
assert measurement(None,c,0) is None
print('Both order orientations, original versus retained coverage, shared metrics and unavailable input passed.')
