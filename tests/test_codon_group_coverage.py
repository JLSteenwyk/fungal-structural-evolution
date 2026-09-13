import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_codon_group_coverage import coverage
class Coverage(unittest.TestCase):
    def test_gap_and_ambiguity_are_not_called(self):
        r=coverage(['ATGAAA','ATGAAG','ATG---','ATGAAN','ATGAAA'])
        self.assertEqual(r['fully_called_columns'],1)
        self.assertEqual(r['columns_at_least_80pct_called'],1)
    def test_variable_covered_column(self):
        r=coverage(['AAA','AAA','AAG','AAG','---'])
        self.assertEqual(r['variable_columns_at_least_80pct_called'],1)
    def test_empty_group(self):self.assertEqual(coverage([])['taxa'],0)
if __name__=='__main__':unittest.main()
