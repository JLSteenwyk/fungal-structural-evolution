#!/usr/bin/env python3
"""Known path/membership and corrupted-export checks for the independent reconstruction."""
import json
from readback_duplication_sister_references import index_tree, reconstruct, check_row

row = dict(gene_node='n1', taxon_id='F1', gene_a='F1_a', gene_b='F1_b')
models = {'F2_r': ('M2', 1), 'F3_s': ('M3', 1), 'F1_r': ('M1', 1)}

def check(tree, dups=(), lookup=models):
    return reconstruct(row, index_tree(tree), set(dups), lookup)

tree = '((F1_a:1,F1_b:2)n1:3,F2_r:4)n0;'
r = check(tree)
assert r['status'] == 'provisional_reference_available'
assert [r[k] for k in ['reference_distance_from_duplicate_node','distance_a_to_reference',
                       'distance_b_to_reference','duplicate_pair_sequence_distance']] == [7,8,9,3]
assert r['nearest_reference_genes'] == ['F2_r']
assert check(tree, ['n0'])['status'] == 'parent_reported_duplication'
assert check(tree, lookup={})['status'] == 'no_modeled_nonfocal_sister'
assert check(tree.replace('F2_r','F1_r'))['status'] == 'focal_taxon_in_sister_clade'
assert check('(F1_a:1,F1_b:2)n1;')['status'] == 'no_parent'
assert check('((F1_a:1,F1_b:2)n1:3,F2_r:4,F3_s:4)n0;')['status'] == 'parent_not_bifurcating'
tied = check('((F1_a:1,F1_b:2)n1:3,(F2_r:1,F3_s:1)n2:3)n0;')
assert tied['nearest_reference_genes'] == ['F2_r','F3_s']
assert tied['chosen_reference_gene'] == 'F2_r'
# Deep ancestors must not enter the immediate-parent sister set or distances.
nested = check('(((F1_a:1,F1_b:2)n1:3,F2_r:4)n0:1000000000000000,F3_s:1)n9;')
assert nested == r
observed = {k: json.dumps(v) if isinstance(v,list) else str(v) for k,v in r.items()}
assert check_row(observed, r) == 0
for field, bad in [('nearest_reference_genes','["F3_s"]'),('reference_model','wrong'),
                   ('sister_genes','2'),('parent_node','n9'),('status','no_parent'),
                   ('distance_a_to_reference','8.01')]:
    damaged = dict(observed, **{field:bad})
    try:
        check_row(damaged, r)
    except ValueError:
        pass
    else:
        raise AssertionError('Accepted corruption: '+field)
for newick in ['((F1_a:1,F1_b:2)n1:-1,F2_r:4)n0;',
               '((F1_a:1,F1_a:2)n1:3,F2_r:4)n0;']:
    try:
        check(newick)
    except ValueError:
        pass
    else:
        raise AssertionError('Accepted malformed tree')
print('Passed known distances, ties, six statuses, large-ancestor isolation, six corruptions and malformed trees.')
