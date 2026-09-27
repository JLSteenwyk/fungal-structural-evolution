#!/usr/bin/env python3
"""Check explicit guide eligibility and model-pair identity compatibility."""
from prepare_background_measurements import eligibility,pair_identity
r=dict(profile_candidate_status='cross_taxon_unreported_candidate',mafft_candidate_status='cross_taxon_unreported_candidate',profile_native_ortholog='1',mafft_native_ortholog='1',both_parents_unreported='1')
assert all(eligibility(r).values())
r['mafft_native_ortholog']='0';x=eligibility(r);assert x['candidate_and_native_ortholog_either']==1 and x['candidate_and_native_ortholog_both']==0
r['profile_candidate_status']='reported_duplication';assert eligibility(r)['candidate_and_native_ortholog_either']==0
assert pair_identity(('M',1),('N',2))==pair_identity(('N',2),('M',1))
assert pair_identity(('M',1),('N',2))[0]!=pair_identity(('M',2),('N',2))[0]
print('Guide/orthology eligibility, missing qualification, pair symmetry and version sensitivity passed.')
