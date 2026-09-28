#!/usr/bin/env python3
"""Check coordinate encoding against independent weighted feature fixtures."""
import json
from pathlib import Path
import unittest
from ancestral_extant_alignment_geometry import project, signature, signature_distance
from check_baliphy_alignment_distance_semantics import distance


class GeometryTests(unittest.TestCase):
    def test_all_independent_metric_fixtures(self):
        receipt = json.loads(Path('metadata/baliphy_alignment_distance_semantics_20260928.json').read_text())
        count = 0
        for group in receipt['evidence']:
            samples = group['alignments']
            observed = {t:s.replace('-','') for t,s in samples[0].items()}
            signatures = [signature(project(s, observed)) for s in samples]
            for i,a in enumerate(samples):
                for j,b in enumerate(samples):
                    self.assertEqual(signature_distance(signatures[i],signatures[j]),distance(a,b))
                    count += 1
        self.assertEqual(count,733)

    def test_ambiguous_input_position_preserved(self):
        observed = {'a':'AX','b':'A'}
        first = project({'a':'AC','b':'A-'},observed)
        second = project({'a':'AG','b':'A-'},observed)
        self.assertEqual(first,second)
        self.assertEqual(signature_distance(signature(first),signature(second)),0)
        moved = project({'a':'AX','b':'-A'},observed)
        self.assertEqual(signature_distance(signature(first),signature(moved)),6)

    def test_internal_only_column_removed(self):
        observed = {'a':'A','b':'A'}
        self.assertEqual(project({'a':'A-','b':'A-','ancestor':'AA'},observed),{'a':'A','b':'A'})

    def test_invalid_observed_inputs(self):
        for sample in [{'a':'G'}, {'a':'AA'}, {'a':'?'}, {'z':'A'}]:
            with self.assertRaises(ValueError):
                project(sample,{'a':'A'})


if __name__ == '__main__':
    unittest.main()
