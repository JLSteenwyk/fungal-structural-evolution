import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from align_tfiib_repeat_pairs import choose_pair

def hit(start,end,coverage=90):return dict(pfam_accession='PF00382.25',alignment_start=start,alignment_end=end,hmm_start=1,hmm_end=coverage,hmm_length=100)
class PairSelection(unittest.TestCase):
    def test_order_is_positional(self):
        pair,status=choose_pair([hit(200,290),hit(10,100)])
        self.assertEqual(pair[0]['alignment_start'],10)
    def test_no_best_two_selection_from_three(self):
        self.assertIsNone(choose_pair([hit(1,90),hit(100,190),hit(200,290)])[0])
    def test_overlap_and_low_coverage(self):
        self.assertIsNone(choose_pair([hit(1,100),hit(100,190)])[0])
        self.assertIsNone(choose_pair([hit(1,90,69),hit(100,190)])[0])
if __name__=='__main__':unittest.main()
