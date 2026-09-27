#!/usr/bin/env python3
"""Cross-check independent full-scenario enumeration on randomized candidate sets."""
import random
from select_background_controls import SCENARIOS,candidate,choose
from readback_background_control_selection import scenarios,reconstruct
assert SCENARIOS==scenarios()
rng=random.Random(92727)
for case in range(100):
    d=0. if case==0 else rng.uniform(.01,3);t=dict(sequence_distance=d);backgrounds={};edges=[]
    for i in range(8):
        bid=str(i);db=0. if d==0 else d*rng.uniform(.5,2)
        backgrounds[bid]=dict(node_id=bid,sequence_distance=db,both_guides=int(i%2==0),both_unreported_parents=int(i%4==0))
        edge=dict(background_id=bid,focal_taxon=str(i%2),shared_genes=str(i==7 and case%2==0),shared_models='0',shared_sequences='0')
        for o in [0,1]:
            edge.update({f'length_ratio_{o}':str(rng.uniform(1,1.6)),f'plddt_difference_{o}':str(rng.uniform(0,16)),f'lowconf_difference_{o}':str(rng.uniform(0,.21))})
        edges.append(edge)
    candidates=[candidate(t,backgrounds[e['background_id']],e) for e in edges]
    for s in SCENARIOS:assert choose(candidates,s,d)==reconstruct(edges,backgrounds,d,s)
# Equal candidates choose lexicographically first background and then order zero.
t=dict(sequence_distance=1);backgrounds={};edges=[]
for bid in ['b','a']:
    backgrounds[bid]=dict(node_id=bid,sequence_distance=1,both_guides=1,both_unreported_parents=1)
    e=dict(background_id=bid,focal_taxon='1',shared_genes='0',shared_models='0',shared_sequences='0')
    for o in [0,1]:e.update({f'length_ratio_{o}':'1',f'plddt_difference_{o}':'0',f'lowconf_difference_{o}':'0'})
    edges.append(e)
for s in SCENARIOS:
    r=reconstruct(edges,backgrounds,1,s);assert r['background_id']=='a' and r['endpoint_order']==0 and r['equal_score_candidates']==2
print('Passed 5400 random scenario reconstructions, zero distances, exclusions and all54 exact-tie cases.')
