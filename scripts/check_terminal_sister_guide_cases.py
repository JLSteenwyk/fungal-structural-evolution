#!/usr/bin/env python3
"""Ensure guide comparison preserves changed labels, disagreement and exclusions."""
from compare_terminal_sister_guides import compare_pair, COPIED, CANDIDATE
base={field:'0' for field in COPIED}
base.update(gene_a='F1_a',gene_b='F2_b',taxon_a='F1',taxon_b='F2',family='OG1',gene_node='n1',candidate_status=CANDIDATE,model_coverage='two_distinct_models',model_a='M1',model_b='M2')
other=dict(base,family='OG2',gene_node='n99')
r=compare_pair(base,other)
assert r['presence']=='both' and r['family_label_agrees']==0 and r['modeled_candidate_in_both_guides']==1
assert r['both_parents_unreported']==1
r=compare_pair(base,dict(other,candidate_status='reported_duplication',parent_reported_duplication='1'))
assert r['candidate_class_agrees']==0 and r['modeled_candidate_in_either_guide']==1 and r['modeled_candidate_in_both_guides']==0 and r['both_parents_unreported']==0
for a,b,status in [(base,None,'profile_only'),(None,base,'mafft_only')]:
    r=compare_pair(a,b);assert r['presence']==status and r['modeled_candidate_in_either_guide']==1
r=compare_pair(dict(base,model_coverage='no_models'),None);assert r['modeled_candidate_in_either_guide']==0
r=compare_pair(dict(base,model_coverage='identical_model'),None);assert r['modeled_candidate_in_either_guide']==1
for field in ['gene_a','model_a','taxon_a','model_coverage']:
    try: compare_pair(base,dict(other,**{field:'changed'}))
    except AssertionError: pass
    else: raise AssertionError('Inconsistent identity accepted')
print('Guide label changes, disagreement, missing guides, coverage states and identity rejection cases passed.')
