"""Rigid-motion controls for partner reconstruction used to audit native 3Di."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_3di_features import coordinate_features


class FeatureGeometry(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(934)
        self.ca = np.cumsum(rng.normal(size=(12, 3)) * 3, axis=0)
        self.n = self.ca + [1., .2, .1]
        self.c = self.ca + [.1, 1., .2]
        self.cb = self.ca + [.2, .1, 1.]
        self.cb[4] = np.nan

    def test_rigid_motion_and_missing_cb(self):
        partner, valid, features = coordinate_features(self.ca, self.n, self.c, self.cb)
        rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        moved = [x @ rotation + [31, -19, 2] for x in [self.ca, self.n, self.c, self.cb]]
        p2, v2, f2 = coordinate_features(*moved)
        np.testing.assert_array_equal(partner, p2)
        np.testing.assert_array_equal(valid, v2)
        np.testing.assert_allclose(features, f2, atol=1e-12)

    def test_endpoints_are_not_observations(self):
        partner, valid, features = coordinate_features(self.ca, self.n, self.c, self.cb)
        self.assertFalse(valid[0]); self.assertFalse(valid[-1])
        self.assertTrue(valid[1:-1].all())
        np.testing.assert_array_equal(partner[[0, -1]], [-1, -1])
        self.assertTrue((features[[0, -1]] == 0).all())

    def test_missing_backbone_requires_review(self):
        self.n[5] = np.nan
        with self.assertRaises(ValueError):
            coordinate_features(self.ca, self.n, self.c, self.cb)


if __name__ == '__main__':
    unittest.main()
