import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from map_pfam_active_sites import project,pattern_status,parse_sites


class ActiveSiteProjectionTests(unittest.TestCase):
    def test_insertions_gaps_and_protein_offset(self):
        mapping=project('ACd-E','xx.xx',10,'ACDE')
        self.assertEqual(mapping,{1:(10,'A'),2:(11,'C'),3:(None,'-'),4:(13,'E')})
        self.assertEqual(pattern_status(mapping,[(1,'A'),(4,'E')]),'complete_pattern_conserved')
        self.assertEqual(pattern_status(mapping,[(3,'H')]),'pattern_has_unaligned_residue')
        self.assertEqual(pattern_status(mapping,[(2,'H')]),'pattern_residues_not_all_conserved')
    def test_changed_source_rejected(self):
        with self.assertRaises(ValueError):project('ACE','xxx',1,'ADE')
        with self.assertRaises(ValueError):pattern_status({1:(1,'A')},[(2,'A')])
    def test_pinned_file_selenocysteine(self):
        p=Path(__file__).resolve().parents[1]/'data/pfam/38.2/active_site.dat'
        patterns=parse_sites(p)
        self.assertEqual(len(patterns),835)
        self.assertTrue(any(aa=='U' for rows in patterns.values() for _,sites in rows for _,aa in sites))


if __name__=='__main__':unittest.main()
