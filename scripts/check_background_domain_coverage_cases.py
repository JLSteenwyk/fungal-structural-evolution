#!/usr/bin/env python3
"""Boundary and exclusion cases for the domain coverage sensitivity screen."""
from screen_background_domain_coverage import reasons

ok=['within_printed_rounding']*2
assert reasons(ok,[70,70],[100,100],[100,100],30,.7)==[]
assert reasons(ok,[70,69],[100,100],[100,100],30,.7)==['low_original_interval_coverage']
assert reasons(ok,[30,30],[100,100],[30,30],30,.7)==['low_original_interval_coverage']
assert 'numeric_discrepancy' in reasons([ok[0],'outside_printed_rounding'],[70,70],[100,100],[100,100],30,.7)
assert reasons(['input_unavailable']*2,[],[30,30],[0,0],30,.7)==['input_unavailable']
assert 'short_alignment' in reasons(ok,[2,2],[30,30],[3,3],30,.7)
for ns in [[30],[31,31]]:
    try:reasons(ok,ns,[30,30],[30,30],30,.7)
    except ValueError:pass
    else:raise AssertionError('Missing order or impossible count accepted')
# Reversing orders and endpoints cannot change the screen decision.
assert reasons(ok,[35,39],[50,60],[40,50],30,.7)==reasons(ok,[39,35],[60,50],[50,40],30,.7)
print('Coverage boundaries, original denominators, quarantine, order symmetry and impossible/missing count rejection passed.')
