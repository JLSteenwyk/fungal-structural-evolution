import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from retrieve_assembly_statistics import parse_statistics
from measure_external_assemblies import n50


class AssemblyStatistics(unittest.TestCase):
    def test_record_n50_boundary_and_invalid_lengths(self):
        self.assertEqual(n50([5, 3, 2]), (5, 1))
        self.assertEqual(n50([4, 3, 3]), (3, 2))
        for lengths in [[], [0, 3], [-1, 3]]:
            with self.assertRaises(ValueError):
                n50(lengths)

    def test_separates_scopes_and_preserves_missingness(self):
        text = ('# GenBank assembly accession: GCA_123456789.1\n'
                'all\tall\tall\tall\ttotal-length\t1000\n'
                'all\tall\tall\tall\ttotal-gap-length\t0\n'
                'all\tall\tall\tall\tcontig-N50\t100\n'
                'Primary Assembly\tall\tall\tall\ttotal-length\t900\n'
                'Primary Assembly\tall\tall\tall\tscaffold-N50\t300\n')
        r = parse_statistics(text, 'GCA_123456789.1')
        self.assertEqual(r['all_assembly_metrics']['total-gap-length'], 0)
        self.assertNotIn('total-gap-length', r['primary_assembly_metrics'])
        self.assertNotIn('contig-N50', r['primary_assembly_metrics'])
        self.assertEqual(r['primary_assembly_metrics']['total-length'], 900)
        for bad in [text.replace('contig-N50\t100', 'contig-N50\t1001'),
                    text + 'all\tall\tall\tall\ttotal-length\t1000\n',
                    text.replace('total-length\t1000', 'total-length\tnan')]:
            with self.assertRaises(ValueError):
                parse_statistics(bad, 'GCA_123456789.1')
        with self.assertRaises(ValueError):
            parse_statistics(text, 'GCA_123456789.2')


if __name__ == '__main__':
    unittest.main()
