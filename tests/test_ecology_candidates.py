import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from import_ecology_candidates import match_genus


class EcologyCandidates(unittest.TestCase):
    def test_genus_evidence_does_not_resolve_ambiguous_names_or_outgroups(self):
        genera = {'Example': [{'primary_lifestyle': 'plant_pathogen'}], 'Duplicate': [{}, {}]}
        status, matches = match_genus({'study_role': 'ingroup', 'species_name': 'Example species'}, genera)
        self.assertEqual(status, 'genus_candidate_requires_species_verification')
        self.assertEqual(len(matches), 1)
        self.assertEqual(match_genus({'study_role': 'ingroup', 'species_name': 'Duplicate species'}, genera)[0], 'ambiguous_genus_rows')
        self.assertEqual(match_genus({'study_role': 'ingroup', 'species_name': 'Renamed species'}, genera)[0], 'no_exact_genus_match')
        self.assertEqual(match_genus({'study_role': 'outgroup', 'species_name': 'Example species'}, genera), ('outgroup_not_assigned', []))


if __name__ == '__main__':
    unittest.main()
