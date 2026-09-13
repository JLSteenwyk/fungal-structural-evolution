import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from validate_published_outgroup_cds import translate_match

class CompleteTranslation(unittest.TestCase):
    def test_terminal_stop_only_is_allowed(self):
        self.assertEqual(translate_match('ATGAAATAA','MK'),('exact_translation',True))
        self.assertEqual(translate_match('ATGAAATTT','MKF'),('exact_translation',False))
        self.assertEqual(translate_match('ATGTAAAAA','MK')[0],'translation_mismatch')

    def test_partial_codon_is_not_trimmed_to_force_agreement(self):
        self.assertEqual(translate_match('ATGAAAT','MK')[0],'non_triplet_cds_length')

    def test_second_terminal_stop_is_not_silently_removed(self):
        self.assertEqual(translate_match('ATGTAATAA','M')[0],'translation_mismatch')

if __name__=='__main__':unittest.main()
