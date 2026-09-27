#!/usr/bin/env python3
"""Direct descendant-pair fixture checks, including discrepancies and unary nodes."""
from validate_duplication_structure_candidates import index_tree,inspect_node
n=index_tree('((a:1,b:2)n1:3,c:4)n0;')
assert inspect_node(n,'n1',['b','a'])==('exact_reported_pair',True)
assert inspect_node(n,'n0',['a','b'])==('descendant_pair_mismatch',False)
assert inspect_node(n,'absent',['a','b'])==('missing_reported_node',False)
n=index_tree('((a:1)n2:1,b:2)n1;')
assert inspect_node(n,'n1',['a','b'])==('exact_reported_pair',False)
for text in ['(a:1,a:2)n0;','(a:1,b:2);']:
    try:index_tree(text)
    except ValueError:pass
    else:raise AssertionError('Ambiguous labels accepted')
print('Six tree-node checks passed')
