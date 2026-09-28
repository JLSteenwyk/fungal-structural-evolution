#!/usr/bin/env python3
import unittest
import numpy as np
from readback_fastml_variant_roundoff import boundary_error, BOUNDARY_TOLERANCE


class BoundaryTests(unittest.TestCase):
    def test_probabilities_in_range(self):
        for v in [0., 1., .5, np.nextafter(1., 0.)]:
            self.assertEqual(boundary_error(v), 0.)

    def test_excursions_are_reported_not_clipped(self):
        for v in [np.nextafter(1., 2.), 1. + BOUNDARY_TOLERANCE, -BOUNDARY_TOLERANCE]:
            self.assertEqual(boundary_error(v), max(-v, v - 1.))
            self.assertGreater(boundary_error(v), 0.)

    def test_beyond_tolerance_is_rejected(self):
        for v in [1. + 2 * BOUNDARY_TOLERANCE, -2 * BOUNDARY_TOLERANCE, 1.000001, -.001]:
            with self.assertRaises(ValueError):
                boundary_error(v)

    def test_nonfinite_is_rejected(self):
        for v in [float('nan'), float('inf'), -float('inf')]:
            with self.assertRaises(ValueError):
                boundary_error(v)


if __name__ == '__main__':
    unittest.main(verbosity=2)
