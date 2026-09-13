import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


mapping = load('map_marker_structures')
retrieval = load('retrieve_matched_models')


class StructuralIntegration(unittest.TestCase):
    def test_insertions_advance_residue_coordinates_but_gaps_do_not(self):
        self.assertEqual(mapping.residue_positions('aaA.-BC'), [1, 2, 3, None, None, 4, 5])

    def test_marker_priority_spreads_taxa_and_preserves_full_queue(self):
        links = [{'taxon_id': t, 'protein_id': p, 'sequence_sha256': p, 'uniprot_accession': a}
                 for t, p, a in [('A', 'p1', 'A1'), ('A', 'p2', 'A2'), ('B', 'p3', 'B1'),
                                  ('B', 'p4', 'B2'), ('A', 'bulk', 'C1'), ('B', 'bulk', 'C1')]]
        keys = {(r['taxon_id'], r['protein_id'], r['sequence_sha256']) for r in links if r['protein_id'] != 'bulk'}
        ordered, count, taxa = retrieval.prioritize(links, keys, set())
        self.assertEqual(ordered, ['A1', 'B1', 'A2', 'B2', 'C1'])
        self.assertEqual((count, taxa), (4, 2))
        remaining, _, _ = retrieval.prioritize(list(reversed(links)), keys, {'A1'})
        self.assertEqual(set(remaining), {'A2', 'B1', 'B2', 'C1'})
        self.assertEqual(len(remaining), len(set(remaining)))


if __name__ == '__main__':
    unittest.main()
