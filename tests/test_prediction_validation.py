"""Detect prediction corruption and confidence scale errors at the artifact boundary."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_marker_predictions import validate_prediction


class PredictionValidation(unittest.TestCase):
    def setUp(self):
        self.pdb = ''.join(f'ATOM  {i:5d}  CA  {aa:3s} A{i:4d}    '
                           f'{i * 3.8:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{p:6.2f}\n'
                           for i, aa, p in [(1, 'ALA', 85.), (2, 'GLY', 70.)])
        self.pae = np.array([[.2, 1.], [2., .3]])
        self.plddt = np.array([85., 70.])

    def test_valid_and_directional_pae(self):
        validate_prediction('AG', self.pdb, self.pae, self.plddt, 31.)

    def test_sequence_mismatch(self):
        with self.assertRaises(ValueError):
            validate_prediction('AA', self.pdb, self.pae, self.plddt, 31.)

    def test_unscaled_confidence(self):
        with self.assertRaises(ValueError):
            validate_prediction('AG', self.pdb, self.pae, self.plddt / 100, 31.)

    def test_truncated_pae(self):
        with self.assertRaises(ValueError):
            validate_prediction('AG', self.pdb, self.pae[:1], self.plddt, 31.)

    def test_nonfinite_pae(self):
        self.pae[0, 1] = np.nan
        with self.assertRaises(ValueError):
            validate_prediction('AG', self.pdb, self.pae, self.plddt, 31.)

    def test_invalid_numbering(self):
        with self.assertRaises(ValueError):
            validate_prediction('AG', self.pdb.replace('A   2', 'A   3'), self.pae, self.plddt, 31.)


if __name__ == '__main__':
    unittest.main()
