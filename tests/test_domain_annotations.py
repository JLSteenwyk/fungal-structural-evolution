import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from summarize_marker_domains import parse_hit, overlap_size


class DomainAnnotations(unittest.TestCase):
    def test_hmmsearch_orientation_and_inclusive_coordinates(self):
        metadata = {'PF00001.1': {'ML': '80', 'ID': 'Example', 'TP': 'Domain',
                                'sequence_ga': 20., 'domain_ga': 15.}}
        line = 'Sabc - 100 Example PF00001.1 80 1e-20 45 0 1 1 1e-19 1e-18 42 0 2 70 10 78 8 82 0.98 Example domain'
        r = parse_hit(line, metadata, {'Sabc': 100})
        self.assertEqual((r['alignment_start'], r['alignment_end']), (10, 78))
        self.assertEqual(r['hmm_coverage'], 69 / 80)
        self.assertEqual(r['pfam_type'], 'Domain')
        for bad in [line.replace('8 82', '11 82'), line.replace('42 0 2', '12 0 2'),
                    line.replace('0.98', 'nan'), line.replace('Sabc', 'Example')]:
            with self.assertRaises(ValueError):
                parse_hit(bad, metadata, {'Sabc': 100})

    def test_overlapping_domains_keep_inclusive_boundary(self):
        self.assertEqual(overlap_size(1, 10, 10, 20), 1)
        self.assertEqual(overlap_size(1, 10, 11, 20), 0)
        self.assertEqual(overlap_size(1, 100, 20, 30), 11)


if __name__ == '__main__':
    unittest.main()
