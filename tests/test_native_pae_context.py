import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from qualify_native_pae import context_maximum

class DirectionalContext(unittest.TestCase):
    def test_reverse_direction_and_terminal_missingness(self):
        valid = np.array([False, True, True, True, True, False])
        partners = np.array([0, 5, 5, 2, 2, 0])
        pae = np.zeros((6, 6))
        pae[5, 0] = 17
        result = context_maximum(valid, partners, pae)
        self.assertEqual(result[1], 17)
        self.assertTrue(np.isnan(result[[0, 5]]).all())
        self.assertEqual(result[2], 0)

    def test_out_of_bounds_partner_rejected(self):
        with self.assertRaises(ValueError):
            context_maximum(np.array([False, True, False]), np.array([0, 0, 0]), np.zeros((3, 3)))

if __name__ == '__main__':
    unittest.main()
