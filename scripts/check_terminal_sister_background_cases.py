#!/usr/bin/env python3
"""Known-tree cases for the terminal sister inventory's eligibility distinctions."""
from inventory_terminal_sister_backgrounds import pairs
models = {'F1_a': ('M1', 1), 'F2_b': ('M2', 1), 'F1_c': ('M1', 1)}
taxa = {'F1', 'F2', 'F3'}
def run(tree, reported=()):
    return list(pairs(tree, set(reported), models, taxa))
r = run('((F2_b:2,F1_a:1)n1:3,F3_c:4)n0;')[0]
assert (r['gene_a'], r['gene_b'], r['sequence_pair_distance']) == ('F1_a', 'F2_b', 3)
assert r['candidate_status'] == 'cross_taxon_unreported_candidate' and r['model_coverage'] == 'two_distinct_models'
assert run('(F1_a:0,F2_b:0)n0;', ['n0'])[0]['candidate_status'] == 'reported_duplication'
r = run('((F1_a:1,F2_b:2)n1:3,F3_c:4)n0;', ['n0'])[0]
assert r['parent_reported_duplication'] == 1 and r['node_reported_duplication'] == 0
r = run('(F1_a:1,F1_c:2)n0;')[0]
assert r['candidate_status'] == 'same_taxon_unreported' and r['model_coverage'] == 'identical_model'
assert run('(F1_a:1,F3_x:2)n0;')[0]['model_coverage'] == 'one_model'
assert run('(F1_x:1,F3_x:2)n0;')[0]['model_coverage'] == 'no_models'
assert run('(F1_a:1,F2_b:2,F3_x:3)n0;') == []
assert len(run('((F1_a:1,F2_b:2)n1:3,(F1_x:1,F3_x:2)n2:4)n0;')) == 2
for tree in ['(F1_a:-1,F2_b:2)n0;', '(F1_a:1,F1_a:2)n0;', '(F9_x:1,F2_b:2)n0;']:
    try:
        run(tree)
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid tree accepted')
print('Passed ordering, distances, duplication flags, coverage states, same-taxon, polytomy, multiple-cherry and invalid-tree cases.')
