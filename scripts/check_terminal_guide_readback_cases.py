#!/usr/bin/env python3
"""Exercise independent sorted-cursor merge and compare known candidate flags."""
import json
from readback_terminal_sister_guides import merge_sources,expected_row
fields=['family','gene_node','candidate_status','model_coverage','parent_reported_duplication','sequence_tip_a','sequence_tip_b','sequence_pair_distance','model_a','version_a','model_b','version_b']
r={k:'0' for k in fields};r.update(gene_a='A',gene_b='B',taxon_a='T1',taxon_b='T2',candidate_status='cross_taxon_unreported_candidate',model_coverage='identical_model')
s=dict(r,family='changed',candidate_status='reported_duplication')
out=expected_row(r,s);assert out['candidate_class_agrees']=='0' and out['family_label_agrees']=='0' and out['modeled_candidate_in_either_guide']=='1' and out['modeled_candidate_in_both_guides']=='0'
assert expected_row(None,r)['presence']=='mafft_only'
a=[('A','B',json.dumps(r)),('C','D',json.dumps(dict(r,gene_a='C',gene_b='D')))]
b=[('A','B',json.dumps(s)),('E','F',json.dumps(dict(r,gene_a='E',gene_b='F')))]
joined=list(merge_sources(a,b));assert len(joined)==3 and joined[0]==(r,s) and joined[1][1] is None and joined[2][0] is None
assert list(merge_sources([],[]))==[]
assert len(list(merge_sources(a,[])))==2
print('Independent merge intersection/left/right/empty cases and disagreement/identical-model flags passed.')
