#!/usr/bin/env python3
"""Synthetic checks for root-HOG parent-clade tracing."""
from trace_root_hog_omissions import trace_family


def main():
    tree = '((a,b)n1,(c,d)n2,e)n0;'
    groups = {'n1': {'a'}, 'n2': {'c', 'd'}}
    row, = trace_family(tree, groups, {'b'}, {'e'})
    assert row == dict(gene_label='e', inside_emitted_parent_clade=False,
                      enclosing_hog_parent='', immediate_gene_tree_parent='n0',
                      immediate_parent_descendant_genes=5,
                      first_ancestor_with_hog_descendants='n0',
                      edges_to_first_hog_join=1, hogs_below_first_join=2)
    bad = [
        (tree, {'n1': {'a', 'c'}}, {'b'}, {'e'}),
        (tree, {'n0': {'a', 'b', 'c', 'd', 'e'}, 'n1': {'a', 'b'}}, set(), {'e'}),
        (tree, groups, {'b'}, {'absent'}),
        (tree, groups, {'b'}, {'a'}),
    ]
    for args in bad:
        try:
            trace_family(*args)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid input accepted')
    row, = trace_family('((a,b)n1,(e,f)n2)n0;', {'n1': {'a', 'b'}}, set(), {'e'})
    assert row['edges_to_first_hog_join'] == 2
    assert row['first_ancestor_with_hog_descendants'] == 'n0'
    print('Passed exact placement, multi-edge join and four rejection checks')


if __name__ == '__main__':
    main()
