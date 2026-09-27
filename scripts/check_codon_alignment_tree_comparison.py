#!/usr/bin/env python3
"""Check split orientation, topology changes, path sums and zero-length edges."""
from io import StringIO
import math

from Bio import Phylo
from compare_codon_alignment_trees import tree_map, pair_distances


def parse(text):
    return Phylo.read(StringIO(text), 'newick')


def main():
    first = parse('(A:1,B:2,(C:3,D:4)80/90:5);')
    reoriented = parse('(C:3,D:4,(B:2,A:1)80/90:5);')
    changed = parse('(A:1,C:3,(B:2,D:4)70/85:5);')
    zero = parse('(A:1,B:2,(C:3,D:4)80/90:0);')
    taxa, edges = tree_map(first)
    assert tree_map(reoriented) == (taxa, edges)
    _, changed_edges = tree_map(changed)
    internal = lambda e: {s for s in e if len(s) > 1}
    assert len(internal(edges) ^ internal(changed_edges)) == 2
    assert internal(tree_map(zero)[1]) == internal(edges)
    assert pair_distances(taxa, edges)[('A', 'D')] == 10
    assert pair_distances(taxa, edges)[('C', 'D')] == 7
    for tree in [first, reoriented, changed, zero]:
        names, branches = tree_map(tree)
        for (a, b), value in pair_distances(names, branches).items():
            assert math.isclose(value, tree.distance(a, b), abs_tol=1e-12)
    for invalid in ['((A:1,B:2)80/90:2,(C:3,D:4)80/90:3);',
                    '(A:1,A:2,(C:3,D:4)80/90:5);',
                    '(A:1,B:2,(C:3,D:4)80/90:-1);']:
        try:
            tree_map(parse(invalid))
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid tree accepted')
    print('Passed root orientation, changed quartet, zero-edge, direct-path and invalid-grid checks')


if __name__ == '__main__':
    main()
