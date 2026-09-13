import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_genus_tree_splits import incompatible

class SplitTests(unittest.TestCase):
    def test_nested_and_complementary_splits_compatible(self):
        u=set('ABCDEF');g=set('ABC')
        for a in [set('AB'),g,u-g,set('DE')]:self.assertFalse(incompatible(a,g,u))
    def test_four_intersections_required_and_root_independent(self):
        u=set('ABCDEF');g=set('ABC');a=set('AD')
        for side in [a,u-a]:
            for group in [g,u-g]:self.assertTrue(incompatible(side,group,u))

if __name__=='__main__':unittest.main()
