import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_3di_feature_overlap import shared_feature_edges, component_sizes


class FeatureOverlapTests(unittest.TestCase):
    def test_remote_shared_residue_is_linked_without_double_count(self):
        features = {1: {1, 2, 3, 90, 91, 92}, 80: {90, 91, 92, 100, 101, 102}, 20: {20, 21, 22}}
        edges = shared_feature_edges(features)
        self.assertEqual(edges, {(1, 80)})
        self.assertEqual(component_sizes(features, edges), [2, 1])

    def test_connected_chain_is_not_count_of_independent_observations(self):
        features = {1: {1, 2}, 2: {2, 3}, 3: {3, 4}}
        self.assertEqual(component_sizes(features, shared_feature_edges(features)), [3])


if __name__ == '__main__':
    unittest.main()
