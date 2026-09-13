import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from map_experimental_ca_residues import classify_atoms


class CAMappingTests(unittest.TestCase):
    def test_missing_alternate_and_modified(self):
        a={'label_comp_id':'ALA','label_alt_id':'.','Cartn_x':'1','Cartn_y':'2','Cartn_z':'3','occupancy':'1'}
        self.assertEqual(classify_atoms([a],'A'),'unambiguous_full_occupancy_CA')
        self.assertEqual(classify_atoms([],'A'),'no_CA_observation')
        self.assertEqual(classify_atoms([a,a],'A'),'multiple_CA_records')
        self.assertEqual(classify_atoms([a|{'occupancy':'.5'}],'A'),'alternate_or_partial_occupancy')
        self.assertEqual(classify_atoms([a|{'label_comp_id':'MSE'}],'M'),'nonstandard_or_mismatching_monomer')
    def test_nonfinite_and_zero_occupancy(self):
        a={'label_comp_id':'ALA','label_alt_id':'.','Cartn_x':'nan','Cartn_y':'2','Cartn_z':'3','occupancy':'1'}
        self.assertEqual(classify_atoms([a],'A'),'invalid_coordinates_or_occupancy')
        self.assertEqual(classify_atoms([a|{'Cartn_x':'1','occupancy':'0'}],'A'),'invalid_coordinates_or_occupancy')


if __name__=='__main__':unittest.main()
