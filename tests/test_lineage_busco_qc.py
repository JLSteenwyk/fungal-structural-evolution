import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_lineage_busco_qc import validate_summary


class LineageBuscoSummaryTests(unittest.TestCase):
    def fixture(self):
        return {'lineage_dataset':{'name':'fungi_odb12.2','creation_date':'2026-05-13','number_of_buscos':'100'},'results':{'n_markers':100,'Single copy BUSCOs':80,'Multi copy BUSCOs':5,'Complete BUSCOs':85,'Fragmented BUSCOs':5,'Missing BUSCOs':10}}
    def check(self,r):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'summary.json';p.write_text(json.dumps(r))
            return validate_summary(p,{'dataset':'fungi_odb12.2','config':{'creation_date':'2026-05-13','number_of_BUSCOs':'100'}})
    def test_valid_partition(self):
        self.assertEqual(self.check(self.fixture())['Complete BUSCOs'],85)
    def test_wrong_panel_or_denominator(self):
        r=self.fixture();r['lineage_dataset']['name']='eukaryota_odb12.2'
        with self.assertRaises(ValueError):self.check(r)
        r=self.fixture();r['results']['n_markers']=125
        with self.assertRaises(ValueError):self.check(r)
    def test_inconsistent_partition(self):
        r=self.fixture();r['results']['Missing BUSCOs']=11
        with self.assertRaises(ValueError):self.check(r)
        r=self.fixture();r['results']['Complete BUSCOs']=84
        with self.assertRaises(ValueError):self.check(r)


if __name__=='__main__':unittest.main()
