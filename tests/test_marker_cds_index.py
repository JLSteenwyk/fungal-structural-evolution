import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from index_marker_cds_sources import unique_record
class Identity(unittest.TestCase):
    def test_identical_multiple_records_are_still_ambiguous(self):
        self.assertEqual(unique_record(['ATG','ATG']),(None,'multiple_source_cds_records'))
    def test_missing_and_unique(self):
        self.assertEqual(unique_record([]),(None,'no_source_cds'))
        self.assertEqual(unique_record(['ATG']),('ATG','unique_source_cds'))
if __name__=='__main__':unittest.main()
