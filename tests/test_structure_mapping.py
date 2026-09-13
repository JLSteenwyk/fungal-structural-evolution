import importlib.util
import unittest
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


mapping = load('map_marker_structures')
retrieval = load('retrieve_matched_models')
comparison = load('compare_marker_structures')


class StructuralIntegration(unittest.TestCase):
    def test_geometry_is_invariant_to_rotation_and_translation(self):
        x = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        y = x @ rotation + [5., 7., -3.]
        result = comparison.geometry(x, y, [1, 4, 7, 10], [2, 5, 8, 11])
        self.assertLess(result['ca_superposition_rmsd_angstrom'], 1e-10)
        self.assertLess(result['local_distance_rms_change_angstrom'], 1e-10)
        mirror = x.copy()
        mirror[:, 0] *= -1
        self.assertGreater(comparison.geometry(x, mirror, [1, 4, 7, 10], [1, 4, 7, 10])['ca_superposition_rmsd_angstrom'], .1)

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
