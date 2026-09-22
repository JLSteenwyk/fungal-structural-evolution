#!/usr/bin/env python3
"""Small adversarial fixtures for unrooted PMSF split readback."""
import io
from Bio import Phylo
from audit_species_pmsf import tree_edges


def parse(text):
    return tree_edges(Phylo.read(io.StringIO(text), 'newick'), ['A', 'B', 'C', 'D'])


def main():
    assert parse('(A:1,B:2,(C:3,D:4):5);') == parse('(D:4,C:3,(B:2,A:1):5);')
    assert set(parse('(A:1,B:2,(C:3,D:4):5);')) != set(parse('(A:1,C:3,(B:2,D:4):5);'))
    for text in ['(A:1,A:2,(C:3,D:4):5);', '(A:1,B:2,(C:3,E:4):5);',
                 '(A:1,B:2,(C:3,D:4));', '(A:1,B:2,(C:3,D:4):-5);',
                 '((A:1,B:2):2,(C:3,D:4):3);']:
        try:
            parse(text)
        except ValueError:
            continue
        raise AssertionError('Accepted malformed or duplicate-edge tree: ' + text)
    print('Passed root-placement invariance, topology discrimination and five malformed-tree checks')


if __name__ == '__main__':
    main()
