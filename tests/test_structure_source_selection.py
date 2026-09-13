import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from map_marker_structures import source_candidates

class SourceSelectionTests(unittest.TestCase):
    def test_higher_confidence_other_pipeline_cannot_replace_primary(self):
        rows = [{'provider': 'GDM', 'tool': 'primary', 'mean_ca_plddt': 70},
                {'provider': 'ATBC', 'tool': 'alternative', 'mean_ca_plddt': 99},
                {'provider': 'GDM', 'tool': 'different-version', 'mean_ca_plddt': 98}]
        self.assertEqual(source_candidates(rows, 'GDM', 'primary'), rows[:1])
        self.assertEqual(source_candidates(rows[1:], 'GDM', 'primary'), [])

if __name__ == '__main__':
    unittest.main()
