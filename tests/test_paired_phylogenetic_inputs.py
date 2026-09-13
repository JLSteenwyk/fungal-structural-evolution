import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from prepare_paired_phylogenetic_inputs import paired_row, site_counts
from prepare_3di_models import validate_model


class PairedInputTests(unittest.TestCase):
    def encoding(self):
        return dict(sequence='ACDE', states='DDDD', valid=np.array([False, True, True, True]),
                    ca_plddt=np.array([90.] * 4), feature_min_plddt=np.array([90., 90., 60., 90.]),
                    feature_max_pae=np.array([0., 5., 5., 11.]))

    def test_coil_is_observed_only_when_valid_and_confident(self):
        a, s, reasons = paired_row('ACDE', [1, 2, 3, 4], {i: (i, 90.) for i in range(1, 5)}, self.encoding())
        self.assertEqual(a, '?C??')
        self.assertEqual(s, '?D??')
        self.assertEqual(reasons['observed'], 1)
        self.assertEqual(sum(reasons.values()), 4)

    def test_wrong_physical_residue_rejected(self):
        with self.assertRaises(ValueError):
            paired_row('ACDE', [2], {2: (3, 90.)}, self.encoding())

    def test_nonfinite_confidence_unobserved(self):
        e = self.encoding()
        e['feature_max_pae'][1] = np.nan
        self.assertEqual(paired_row('ACDE', [2], {2: (2, 90.)}, e)[:2], ('?', '?'))

    def test_missing_values_do_not_add_states(self):
        self.assertEqual(site_counts(['AA?C', 'AA?C', 'AC?D', 'AC??']), (2, 1))

    def test_model_rejects_bad_dimensions_and_frequency_sum(self):
        with self.assertRaises(ValueError):
            validate_model(b'1 2\n')
        text = '\n'.join(' '.join(['1'] * n) for n in range(1, 20)) + '\n' + ' '.join(['.1'] * 20)
        with self.assertRaises(ValueError):
            validate_model(text.encode())


if __name__ == '__main__':
    unittest.main()
