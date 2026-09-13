import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from build_marker_codon_alignments import project_codons
class CodonColumns(unittest.TestCase):
    def test_gap_and_mask_preserve_full_protein_positions(self):
        self.assertEqual(project_codons('ATGAAAGGGTAA','M-KG',1,[2,4]),('---GGG',True))
    def test_nonstandard_code(self):
        self.assertEqual(project_codons('CTG','S',12,[1]),('CTG',False))
    def test_mismatch_outside_selected_columns_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'full_protein_translation_mismatch'):project_codons('ATGAAA','MG',1,[1])
    def test_no_second_stop_or_partial_codon_repair(self):
        for dna in ['ATGTAATAA','ATGA']:
            with self.assertRaises(ValueError):project_codons(dna,'M',1,[1])
    def test_invalid_mask(self):
        with self.assertRaisesRegex(ValueError,'invalid_column_mask'):project_codons('ATG','M',1,[1,1])
if __name__=='__main__':unittest.main()
