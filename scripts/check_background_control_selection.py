#!/usr/bin/env python3
"""Check deterministic choices, qualification/caliper exclusions, ties and exact zero."""
from select_background_controls import candidate,choose,SCENARIOS
assert len(SCENARIOS)==54 and len({s['scenario_id'] for s in SCENARIOS})==54
edge=dict(focal_taxon='1',shared_genes='0',shared_models='0',shared_sequences='0')
for o in [0,1]:edge.update({f'length_ratio_{o}':'1',f'plddt_difference_{o}':'0',f'lowconf_difference_{o}':'0'})
t=dict(sequence_distance=1.)
b=dict(node_id='b',sequence_distance=1.,both_guides=1,both_unreported_parents=1)
c=candidate(t,b,edge);a=candidate(t,dict(b,node_id='a'),edge)
for s in SCENARIOS:
 r=choose([c,a],s,1.);assert r['background_id']=='a' and r['endpoint_order']==0 and r['score']==0 and r['equal_score_candidates']==2 and r['score_gap_to_second']==0
 assert choose([dict(c,shared=True)],s,1.) is None
 if s['focal_only']:assert choose([dict(c,focal=0)],s,1.) is None
 if s['background_set']!='guide_native_ortholog':assert choose([candidate(t,dict(b,both_guides=0,both_unreported_parents=0),edge)],s,1.) is None
 assert choose([],s,1.) is None
z=candidate(dict(sequence_distance=0),dict(b,sequence_distance=0),edge);assert choose([z],SCENARIOS[0],0)['score']==0
try:candidate(dict(sequence_distance=0),b,edge)
except ValueError:pass
else:raise AssertionError('Zero/positive distance accepted')
# Length and confidence fit under different orientations; joint calipers must fail.
e=dict(edge,length_ratio_0='1',plddt_difference_0='40',length_ratio_1='2',plddt_difference_1='0')
bad=candidate(t,b,e);assert all(choose([bad],s,1.) is None for s in SCENARIOS)
# Equal multiplicative bounds preserve original comparisons at the threshold.
boundary=candidate(t,dict(b,sequence_distance=1.5),edge)
for s in SCENARIOS:assert (choose([boundary],s,1.) is not None)==(s['distance_factor']>=1.5)
print('Passed all54 scenarios, deterministic ties, shared/focal/qualification exclusions, empty pools, exact zero, joint-orientation and distance-boundary cases.')
