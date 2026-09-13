"""Mechanical controls for separating global placement from within-domain geometry."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_marker_domains import fit_residuals


class DomainGeometry(unittest.TestCase):
    def test_rigid_transform_invariant(self):
        x = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])
        rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        self.assertLess(fit_residuals(x, x @ rotation + [8, 9, 2]).max(), 1e-20)

    def test_independent_domain_translation(self):
        domain = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])
        x = np.vstack([domain, domain + [20, 0, 0]])
        y = np.vstack([domain, domain + [20, 9, 0]])
        self.assertGreater(fit_residuals(x, y).mean(), 1.)
        self.assertLess(fit_residuals(x[4:], y[4:]).max(), 1e-20)

    def test_reflection_not_accepted_as_rotation(self):
        x = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])
        self.assertGreater(fit_residuals(x, x * [-1, 1, 1]).mean(), .1)


if __name__ == '__main__':
    unittest.main()
