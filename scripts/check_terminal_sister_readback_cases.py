#!/usr/bin/env python3
"""Check independent reconstruction and rejection of altered exported fields."""
from readback_terminal_sister_backgrounds import reconstruct, compare
lookup = {'F1_a': ('M1', 1), 'F2_b': ('M2', 1)}
row = reconstruct('((F2_b:2,F1_a:1)n1:3,F3_c:4)n0;', {'n0'}, lookup, {'F1', 'F2', 'F3'})['n1']
assert row['gene_a'] == 'F1_a' and row['sequence_pair_distance'] == 3
assert row['parent_reported_duplication'] == '1' and row['node_reported_duplication'] == '0'
actual = {k: str(v) for k, v in row.items()}
compare(actual, row)
for field, value in [('gene_a', 'F1_wrong'), ('sequence_pair_distance', '4'), ('version_a', '2'), ('parent_reported_duplication', '0'), ('candidate_status', 'reported_duplication')]:
    changed = dict(actual); changed[field] = value
    try:
        compare(changed, row)
    except AssertionError:
        pass
    else:
        raise AssertionError('Altered field accepted: ' + field)
assert reconstruct('(F1_a:1,F2_b:2,F3_c:3)n0;', set(), lookup, {'F1','F2','F3'}) == {}
print('Independent topology/distance fixture and five exported-field corruption cases passed.')
