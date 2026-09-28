#!/usr/bin/env python3
import unittest
from ancestral_residue_anchors import anchored_columns


class AnchorTests(unittest.TestCase):
    def test_shifted_columns_preserve_extant_identity(self):
        a = anchored_columns({'a': 'AC', 'b': 'A-', 'n': 'AC'},
                             {'a': 'AC', 'b': 'A'}, {'source': 'n'})
        b = anchored_columns({'a': '-AC', 'b': '-A-', 'other': 'GAC'},
                             {'a': 'AC', 'b': 'A'}, {'source': 'other'})
        self.assertEqual(a, b[1:])
        self.assertTrue(b[0]['unanchored'])
        self.assertEqual(b[0]['states'], {'source': 'G'})

    def test_homology_changes_not_hidden_by_same_letters(self):
        observed = {'a': 'AA', 'b': 'A'}
        a = anchored_columns({'a': 'AA', 'b': 'A-', 'n': 'AA'}, observed, {'node': 'n'})
        b = anchored_columns({'a': 'AA', 'b': '-A', 'n': 'AA'}, observed, {'node': 'n'})
        self.assertNotEqual(a[0]['anchors'], b[0]['anchors'])
        self.assertIn(('b', 1), a[0]['anchors'])
        self.assertIn(('b', 1), b[1]['anchors'])

    def test_tip_order_and_runtime_names_do_not_define_identity(self):
        a = anchored_columns({'a': 'AC', 'b': 'A-', 'n1': 'AV'},
                             {'a': 'AC', 'b': 'A'}, {'node': 'n1'})
        b = anchored_columns({'b': 'A-', 'n9': 'AV', 'a': 'AC'},
                             {'b': 'A', 'a': 'AC'}, {'node': 'n9'})
        self.assertEqual(a, b)

    def test_unanchored_columns_are_retained_without_invented_identity(self):
        rows = anchored_columns({'a': '--A-', 'n': 'GVAC'}, {'a': 'A'}, {'node': 'n'})
        self.assertEqual([r['unanchored'] for r in rows], [True, True, False, True])
        self.assertEqual([r['anchors'] for r in rows], [(), (), (('a', 1),), ()])
        self.assertEqual(''.join(r['states']['node'] for r in rows), 'GVAC')

    def test_ambiguous_input_allows_sampled_state_but_known_residue_cannot_change(self):
        anchored_columns({'a': 'AC', 'n': 'AA'}, {'a': 'AX'}, {'node': 'n'})
        with self.assertRaises(ValueError):
            anchored_columns({'a': 'AC', 'n': 'AA'}, {'a': 'AA'}, {'node': 'n'})

    def test_invalid_alignment_and_missing_nodes(self):
        for seq in [{'a': 'A', 'n': 'AC'}, {'a': 'A', 'n': '?'}, {'a': 'A'}]:
            with self.assertRaises(ValueError):
                anchored_columns(seq, {'a': 'A'}, {'node': 'n'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
