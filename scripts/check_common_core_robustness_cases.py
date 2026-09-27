#!/usr/bin/env python3
from summarize_common_core_robustness import direction
assert direction([1]*7,8,0)=='incomplete'
assert direction([.2]*8,8,.1)=='a_farther_from_reference'
assert direction([-.2]*8,8,.1)=='b_farther_from_reference'
assert direction([.2]*7+[-.2],8,.1)=='changes_direction_across_alternatives'
assert direction([.1]*8,8,.1)=='touches_or_within_margin_band'
assert direction([.2]*7+[0],8,.1)=='touches_or_within_margin_band'
assert direction([0]*8,8,0)=='touches_or_within_margin_band'
print('Incomplete grids, opposing directions, zero and exact-margin boundaries checked.')
