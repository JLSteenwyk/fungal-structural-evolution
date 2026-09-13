import importlib.util
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

spec = importlib.util.spec_from_file_location('audit', Path(__file__).resolve().parents[1] / 'scripts/verify_sanchytrid_coordinates.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class Coordinates(unittest.TestCase):
    def test_alias_requires_exact_translation(self):
        genomes = {'NODE_1_length_9_cov_1': SeqRecord(Seq('ATGGCTTAA'))}
        aliases = {'NODE_1_length_9': list(genomes)}
        record = SeqRecord(Seq('MA*'), id='p', description='p NODE_1_length_9_cov_0.1:1-9(+)')
        result, cds = audit.verify(record, genomes, aliases)
        self.assertEqual(result['status'], 'exact_translation')
        self.assertEqual(cds, 'ATGGCTTAA')
        record.seq = Seq('MM*')
        self.assertEqual(audit.verify(record, genomes, aliases)[0]['status'], 'translation_mismatch')

    def test_boundaries_and_reverse_strand(self):
        genomes = {'ctg': SeqRecord(Seq('TTAAGCCAT'))}
        record = SeqRecord(Seq('MA*'), id='p', description='p ctg:9-1(-)')
        self.assertEqual(audit.verify(record, genomes, {})[0]['status'], 'exact_translation')
        record.description = 'p ctg:10-2(-)'
        self.assertEqual(audit.verify(record, genomes, {})[0]['status'], 'coordinates_out_of_bounds')


if __name__ == '__main__':
    unittest.main()
