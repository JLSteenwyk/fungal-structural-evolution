import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from benchmark_tree_paths_geometry import path_mask, path_interval


class PathTests(unittest.TestCase):
    def test_splits_separate_endpoints(self):
        splits = [('A',), ('B',), ('C',), ('D',), ('A', 'B')]
        np.testing.assert_array_equal(path_mask(splits, 'A', 'C'), [True, False, True, False, True])
        np.testing.assert_array_equal(path_mask(splits, 'A', 'B'), [True, True, False, False, False])

    def test_joint_covariance_preserved_in_path_interval(self):
        draws = [[0., 1.], [.25, .75], [.75, .25], [1., 0.]]
        np.testing.assert_allclose(path_interval(draws, np.array([True, True])), [1., 1., 1.])
        marginal_upper_sum = np.quantile(np.array(draws), .975, axis=0).sum()
        self.assertGreater(marginal_upper_sum, 1.)
        with self.assertRaises(ValueError):
            path_interval([[0., -1.]], np.array([True, True]))


if __name__ == '__main__':
    unittest.main()
