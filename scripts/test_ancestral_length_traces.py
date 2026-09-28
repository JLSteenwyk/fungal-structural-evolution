#!/usr/bin/env python3
"""Check sparse iteration matching, node identity and duplicate/missing samples."""
import copy
import unittest
from prepare_ancestral_length_diagnostics import length_rows


class LengthTests(unittest.TestCase):
    def setUp(self):
        self.samples = [dict(iteration=i, level=str(level), source_node='n%d' % level,
                             runtime_node='process-local', ungapped_length=100+i+level)
                        for i in [0, 10, 20] for level in range(4)]

    def test_iteration_and_node_identity(self):
        rows, nodes = length_rows(list(reversed(self.samples)), [0, 10, 20])
        self.assertEqual([r['iter'] for r in rows], [0, 10, 20])
        self.assertEqual(rows[2]['length_level3'], 123)
        self.assertEqual(nodes[3], 'n3')

    def test_missing_and_duplicate(self):
        for samples in [self.samples[:-1], self.samples + [self.samples[0]]]:
            with self.assertRaises(AssertionError):
                length_rows(samples, [0, 10, 20])

    def test_changed_node(self):
        samples = copy.deepcopy(self.samples)
        samples[-1]['source_node'] = 'different-node'
        with self.assertRaises(AssertionError):
            length_rows(samples, [0, 10, 20])

    def test_invalid_length(self):
        for value in [-1, 1.5]:
            samples = copy.deepcopy(self.samples)
            samples[0]['ungapped_length'] = value
            with self.assertRaises(AssertionError):
                length_rows(samples, [0, 10, 20])


if __name__ == '__main__':
    unittest.main(verbosity=2)
