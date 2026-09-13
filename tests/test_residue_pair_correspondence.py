import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from realign_codon_review_groups import residue_pairs
class Correspondence(unittest.TestCase):
    def test_positions_count_residues_not_columns(self):
        pairs,diff=residue_pairs('A-CD','AB-D',[1,2,3,4],{1,2,3},{1,2,3})
        self.assertEqual(pairs,{(1,1),(3,3)});self.assertEqual(diff,0)
    def test_source_codon_ambiguity_excluded(self):
        pairs,diff=residue_pairs('AC','AD',[1,2],{1},{1,2})
        self.assertEqual(pairs,{(1,1)});self.assertEqual(diff,0)
    def test_mask_does_not_reset_residue_positions(self):
        pairs,diff=residue_pairs('AC','AD',[2],{1,2},{1,2})
        self.assertEqual(pairs,{(2,2)});self.assertEqual(diff,1)
if __name__=='__main__':unittest.main()
