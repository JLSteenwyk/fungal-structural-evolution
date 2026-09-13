import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from screen_experimental_sequences import classify


class ExperimentalSequenceTests(unittest.TestCase):
    def test_exact_full_and_unique_fragment(self):
        self.assertEqual(classify('ACDE','ACDE'),('exact_full_sequence',[1]))
        self.assertEqual(classify('CDE','AACDEF'),('exact_fragment_unique',[3]))
    def test_repeated_fragment_is_ambiguous(self):
        self.assertEqual(classify('ACA','ACACA'),('exact_fragment_ambiguous',[1,3]))
    def test_construct_and_variant_are_not_full_matches(self):
        self.assertEqual(classify('HHACDE','ACDE')[0],'target_contained_in_longer_construct')
        self.assertEqual(classify('ACDF','ACDE')[0],'requires_alignment_or_variant_review')
    def test_unknown_and_missing_not_exact(self):
        self.assertEqual(classify('ACX','ACX')[0],'noncanonical_sequence_requires_review')
        self.assertEqual(classify('','ACDE')[0],'missing_sequence')


if __name__=='__main__':unittest.main()
