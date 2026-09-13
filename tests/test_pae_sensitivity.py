import json
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from retrieve_marker_pae import validate_pae
from assess_pae_sensitivity import confidence_mask


class PAESensitivity(unittest.TestCase):
    def test_rejects_wrong_dimensions_nonfinite_and_invalid_range(self):
        def payload(matrix):
            return json.dumps([{'predicted_aligned_error': matrix, 'max_predicted_aligned_error': 31.75}])
        self.assertEqual(validate_pae(payload([[0, 1], [2, 0]]), 2).shape, (2, 2))
        for matrix in [[[0]], [[0, -1], [2, 0]], [[0, float('nan')], [2, 0]], [[0, 40], [2, 0]]]:
            with self.assertRaises(ValueError):
                validate_pae(payload(matrix), 2)

    def test_requires_both_directions_in_both_models_at_mapped_residues(self):
        a, b = np.zeros((5, 5)), np.zeros((6, 6))
        # Alignment maps protein residues 2,5 to 1,4; use asymmetric errors.
        a[4, 1] = 12
        self.assertFalse(confidence_mask(a, b, [2, 5], [1, 4], 10)[0, 1])
        self.assertTrue(confidence_mask(a, b, [2, 5], [1, 4], 15)[0, 1])
        a[4, 1] = 0
        b[0, 3] = 11
        mask = confidence_mask(a, b, [2, 5], [1, 4], 10)
        self.assertFalse(mask[0, 1])
        np.testing.assert_array_equal(mask, mask.T)
        with self.assertRaises(ValueError):
            confidence_mask(a, b, [0, 5], [1, 4], 10)


if __name__ == '__main__':
    unittest.main()
