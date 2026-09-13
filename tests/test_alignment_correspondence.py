import itertools
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_alignment_correspondence import agreement


class AlignmentCorrespondence(unittest.TestCase):
    def test_counts_match_explicit_cross_taxon_edges(self):
        profile = {(t, p): p for t in 'ABC' for p in [1, 2]}
        mafft = {('A', 1): 1, ('B', 1): 1, ('C', 1): 2,
                 ('A', 2): 2, ('B', 2): 2, ('C', 2): 3, ('D', 9): 3}
        common = profile.keys() & mafft.keys()
        def edges(mapping):
            return {frozenset((a, b)) for a, b in itertools.combinations(common, 2)
                    if a[0] != b[0] and mapping[a] == mapping[b]}
        pe, me = edges(profile), edges(mafft)
        result, columns = agreement(profile, mafft)
        self.assertEqual(result['profile_edges_on_common_residues'], len(pe))
        self.assertEqual(result['mafft_edges_on_common_residues'], len(me))
        self.assertEqual(result['shared_edges'], len(pe & me))
        self.assertEqual(result['edge_jaccard'], .25)
        self.assertEqual(sum(r['shared_edges'] for r in columns), result['shared_edges'])
        self.assertEqual(result['mafft_retained_residues'], 7)
        self.assertEqual(result['common_retained_residues'], 6)

    def test_column_labels_do_not_affect_agreement(self):
        a = {('A', 1): 1, ('B', 2): 1, ('C', 3): 2}
        b = {('A', 1): 8, ('B', 2): 8, ('C', 3): 11}
        result, _ = agreement(a, b)
        self.assertEqual(result['edge_jaccard'], 1)
        singleton, _ = agreement({('A', 1): 1}, {('A', 1): 3})
        self.assertEqual(singleton['edge_jaccard'], '')


if __name__ == '__main__':
    unittest.main()
