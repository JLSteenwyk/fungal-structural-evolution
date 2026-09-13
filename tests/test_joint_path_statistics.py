import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from append_paired_path_uncertainty import joint_path_statistics


class JointPathTests(unittest.TestCase):
    def test_anticorrelated_branches_have_constant_path(self):
        result = joint_path_statistics([[0., 10.], [10., 0.]], [[0., 20.], [20., 0.]], np.array([[True, True]]))
        for suffix in ['p025', 'median', 'p975']:
            self.assertEqual(result['aa_' + suffix][0], 10.)
            self.assertEqual(result['3di_af_' + suffix][0], 20.)
        self.assertEqual(result['aa_sd'][0], 0.)
        self.assertEqual(result['paired_sampling_covariance'][0], 0.)

    def test_joint_covariance_retains_draw_pairing(self):
        result = joint_path_statistics([[1.], [2.], [3.]], [[3.], [2.], [1.]], np.array([[True]]))
        self.assertEqual(result['paired_sampling_covariance'][0], -1.)
        self.assertEqual(result['aa_median'][0], 2.)

    def test_rejects_mismatched_or_invalid_draws(self):
        for aa, st in [([[1.]], [[1.]]), ([[1.], [2.]], [[1.]]), ([[1.], [-2.]], [[1.], [2.]]), ([[1.], [float('nan')]], [[1.], [2.]])]:
            with self.assertRaises(ValueError):
                joint_path_statistics(aa, st, np.array([[True]]))


if __name__ == '__main__':
    unittest.main()
