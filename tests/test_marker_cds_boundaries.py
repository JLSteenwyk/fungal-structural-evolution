import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_ncbi_marker_boundaries import describe

def part(start,end,strand='+',phase='0',**attrs):
    return dict(seqid='s',start=start,end=end,strand=strand,phase=phase,attributes=attrs)
class Boundaries(unittest.TestCase):
    def test_minus_strand_end_range_marks_five_prime(self):
        r=describe([part(10,18,'-',end_range=['18','.'])],{})
        self.assertTrue(r['five_prime_partial']);self.assertFalse(r['three_prime_partial'])
    def test_internal_boundary_is_not_terminal(self):
        r=describe([part(1,3,end_range=['3','.']),part(8,10)],{})
        self.assertTrue(r['internal_partial_boundary']);self.assertFalse(r['three_prime_partial'])
    def test_regional_code_and_default(self):
        self.assertEqual(describe([part(1,6)],{'s':{'12'}})['translation_table'],'12')
        self.assertEqual(describe([part(1,6)],{})['code_provenance'],'ncbi_documented_default')
    def test_code_disagreement_is_flagged(self):
        r=describe([part(1,6,transl_table=['12'])],{'s':{'1'}})
        self.assertIn('cds_region_code_disagreement',r['annotation_issues'])
    def test_split_codon_phase(self):
        self.assertEqual(describe([part(1,4),part(8,12,phase='2')],{})['annotation_status'],'consistent_ordered_cds_features')
if __name__=='__main__':unittest.main()
