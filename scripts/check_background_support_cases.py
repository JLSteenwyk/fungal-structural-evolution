#!/usr/bin/env python3
"""Check distance-window counts against enumeration and guide set eligibility."""
from assess_background_matching_support import interval_count,supported_sets,SETS
values=[0,0,.5,.8,1,1.25,1.5,2,4]
for distance in [0,.1,.5,1,2,5]:
    for factor in [1.25,1.5,2.0]:assert interval_count(values,distance,factor)==sum(distance/factor<=v<=distance*factor for v in values)
assert interval_count(values,0,2)==2 and interval_count([],1,2)==0
r=dict(profile_candidate_status='cross_taxon_unreported_candidate',profile_native_ortholog='1',mafft_native_ortholog='1',modeled_candidate_in_both_guides='1',both_parents_unreported='1')
assert supported_sets(r,'profile')==list(SETS)
r['both_parents_unreported']='0';assert supported_sets(r,'profile')==list(SETS[:2])
r['mafft_native_ortholog']='0';assert supported_sets(r,'profile')==list(SETS[:1])
r['profile_native_ortholog']='0';assert supported_sets(r,'profile')==[]
print('18 enumerated windows, duplicate zeros, empty pool and three nested eligibility sets passed.')
