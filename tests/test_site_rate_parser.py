import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from estimate_paired_site_rates import rate_rows

class SiteRateParserTest(unittest.TestCase):
    def test_both_documented_and_actual_headers(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rates'
            for header in ['Site Rate Cat C_Rate','Site Rate Category Categorized_rate']:
                p.write_text('# annotation\n'+header+'\n1 0.25 1 0.10\n2 1.5 4 2.0\n')
                r=rate_rows(p,2);self.assertEqual(r[0]['Category'],'1');self.assertEqual(r[1]['Categorized_rate'],'2.0')
    def test_invalid_values_and_grid(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rates'
            for row in ['2 0.2 1 0.1','1 nan 1 0.1','1 0.2 5 0.1','1 -0.2 1 0.1']:
                p.write_text('Site Rate Cat C_Rate\n'+row+'\n')
                with self.assertRaises(ValueError):rate_rows(p,1)
if __name__=='__main__':unittest.main()
