import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_ncbi_cds_translation import compare,code_choice

class StrictTranslation(unittest.TestCase):
    def test_annotated_code_changes_cug_translation(self):
        self.assertEqual(compare('CTG','S',12)[0],'exact_translation')
        self.assertEqual(compare('CTG','S',1)[0],'translation_mismatch')
    def test_missing_code_is_explicit_assumption(self):
        self.assertEqual(code_choice(set()),(1,'table_1_assumption_no_explicit_code'))
        self.assertIsNone(code_choice({'1','12'})[0])
    def test_no_partial_or_initiator_repair(self):
        self.assertEqual(compare('ATGA','M',1)[0],'non_triplet_length')
        self.assertEqual(compare('GTG','M',1)[0],'translation_mismatch')
    def test_only_one_terminal_stop_removed(self):
        self.assertEqual(compare('ATGTAA','M',1),('exact_translation',True))
        self.assertEqual(compare('ATGTAATAA','M',1)[0],'translation_mismatch')

if __name__=='__main__':unittest.main()
