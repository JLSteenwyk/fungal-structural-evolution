import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_paired_branch_resampling import sampled_columns, draw_seed


class ResamplingTests(unittest.TestCase):
    def test_length_bounds_and_within_block_contiguity(self):
        indices = sampled_columns(91, 30, 17)
        self.assertEqual(len(indices), 91)
        self.assertTrue(np.all((0 <= indices) & (indices < 91)))
        for start in [0, 30, 60]:
            np.testing.assert_array_equal(np.diff(indices[start:start+30]) % 91, np.ones(29, dtype=int))

    def test_shared_indices_preserve_missingness(self):
        a = np.array([list('AC??DEF'), list('A?DE?FG')])
        s = np.array([list('DD??DDD'), list('D?DD?DD')])
        indices = sampled_columns(7, 3, draw_seed('marker', 3, 0))
        np.testing.assert_array_equal(a[:, indices] == '?', s[:, indices] == '?')

    def test_determinism_and_full_block(self):
        seed = draw_seed('m', 5, 0)
        np.testing.assert_array_equal(sampled_columns(5, 5, seed), sampled_columns(5, 5, seed))
        self.assertEqual(set(sampled_columns(5, 5, seed)), set(range(5)))
        self.assertNotEqual(seed, draw_seed('m', 5, 1))
        with self.assertRaises(ValueError):
            sampled_columns(5, 6, seed)


if __name__ == '__main__':
    unittest.main()
