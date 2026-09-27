#!/usr/bin/env python3
"""Known-tree checks for provisional reference eligibility and path lengths."""
from inventory_duplication_sister_references import tree_index,reference_summary
r={'gene_node':'n1','taxon_id':'F1','gene_a':'F1_a','gene_b':'F1_b'}
models={'F2_r':('M2',1),'F3_s':('M3',1),'F1_r':('M1',1)};taxa={'F1','F2','F3'}
def check(tree,dups=set(),lookup=models):return reference_summary(r,tree_index(tree),dups,lookup,taxa)
t='((F1_a:1,F1_b:2)n1:3,F2_r:4)n0;';v=check(t)
assert v['status']=='provisional_reference_available' and v['chosen_reference_gene']=='F2_r'
assert (v['reference_distance_from_duplicate_node'],v['distance_a_to_reference'],v['distance_b_to_reference'],v['duplicate_pair_sequence_distance'])==(7,8,9,3)
assert check(t,{'n0'})['status']=='parent_reported_duplication'
assert check(t,lookup={})['status']=='no_modeled_nonfocal_sister'
assert check(t.replace('F2_r','F1_r'))['status']=='focal_taxon_in_sister_clade'
v=check('((F1_a:1,F1_b:2)n1:3,(F2_r:1,F3_s:1)n2:3)n0;');assert v['nearest_reference_genes']=='["F2_r", "F3_s"]' and v['chosen_reference_gene']=='F2_r'
assert check('((F1_a:1,F1_b:2)n1:3,F2_r:4,F3_s:4)n0;')['status']=='parent_not_bifurcating'
assert check('(F1_a:1,F1_b:2)n1;')['status']=='no_parent'
assert check('((F1_a:0,F1_b:0)n1:0,F2_r:0)n0;')['reference_distance_from_duplicate_node']==0
for tree in ['((F1_a:1,F1_b:2)n1:-1,F2_r:4)n0;','((F1_a:1,F1_a:2)n1:3,F2_r:4)n0;']:
    try:check(tree)
    except ValueError:pass
    else:raise AssertionError('Invalid branch/label accepted')
print('Reference eligibility, tied references, four path distances, root/polytomy/nested-duplication/focal-taxon exclusions and invalid trees checked.')
