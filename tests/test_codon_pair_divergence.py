import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_codon_pair_divergence import encode,differences
class PairDifferences(unittest.TestCase):
    def test_same_amino_acid_and_shared_denominator(self):
        r=differences(encode('AAAAAG---',1),encode('AAGAAGAAN',1))
        self.assertEqual(r['shared_called_codons'],2)
        self.assertEqual(r['different_codons_same_amino_acid'],1)
        self.assertEqual(r['different_amino_acids'],0)
        self.assertEqual(r['third_position_difference_fraction'],.5)
    def test_code_specific_amino_acid_change(self):
        r=differences(encode('CTG',12),encode('TCT',12))
        self.assertEqual(r['different_amino_acids'],0)
        r=differences(encode('CTG',1),encode('TCT',1))
        self.assertEqual(r['different_amino_acids'],1)
    def test_no_shared_data_is_not_zero_divergence(self):
        r=differences(encode('---',1),encode('AAA',1))
        self.assertEqual(r['codon_difference_fraction'],'')
    def test_stop_is_rejected(self):
        with self.assertRaises(ValueError):encode('TAA',1)
if __name__=='__main__':unittest.main()
