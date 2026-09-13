import io,sys,unittest
from pathlib import Path
from Bio import Phylo
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_tfiib_domain_class_split import masks


class SplitMaskTests(unittest.TestCase):
    def test_split_independent_of_root_placement(self):
        index={t:i for i,t in enumerate('abcd')};full=15
        sets=[]
        for newick in ['((a:1,b:1):2,c:1,d:1);','(a:1,b:1,(c:1,d:1):2);']:
            cm=masks(Phylo.read(io.StringIO(newick),'newick'),index)
            sets.append({min(v,full^v) for v in cm.values() if v.bit_count()==2})
        self.assertEqual(sets[0],sets[1]);self.assertEqual(sets[0],{3})
    def test_duplicate_tips_and_bad_lengths_rejected(self):
        for newick in ['(a:1,a:1,b:1);','(a:-1,b:1);']:
            with self.assertRaises(ValueError):masks(Phylo.read(io.StringIO(newick),'newick'),{'a':0,'b':1})


if __name__=='__main__':unittest.main()
