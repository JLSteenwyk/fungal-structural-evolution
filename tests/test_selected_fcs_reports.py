import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('fcs',Path(__file__).resolve().parents[1]/'scripts/retrieve_selected_fcs_reports.py')
fcs=importlib.util.module_from_spec(spec);spec.loader.exec_module(fcs)
HEADER='##[["FCS genome report",2,1],{}]\n#'+'\t'.join(fcs.FIELDS)+'\n'


class FCSReports(unittest.TestCase):
    def test_header_only_is_valid_but_empty_file_is_not(self):
        self.assertEqual(fcs.parse_report(HEADER)[1],[])
        with self.assertRaises(ValueError):fcs.parse_report('')

    def test_one_based_union_and_action_separation(self):
        text=HEADER+'s1\t1\t10\t100\tFIX\tx\t50\tx\n'+'s1\t5\t15\t100\tTRIM\tx\t50\tx\n'+'s1\t20\t20\t100\tREVIEW\tx\t0\tx\n'
        rows=fcs.parse_report(text)[1]
        self.assertEqual(fcs.union_bases(rows,{'FIX','TRIM'}),15)
        self.assertEqual(fcs.union_bases(rows,{'REVIEW'}),1)

    def test_invalid_coordinates_and_duplicate_records_rejected(self):
        with self.assertRaises(ValueError):fcs.parse_report(HEADER+'s1\t0\t10\t100\tFIX\tx\t50\tx\n')
        row='s1\t1\t10\t100\tFIX\tx\t50\tx\n'
        with self.assertRaises(ValueError):fcs.parse_report(HEADER+row+row)


if __name__=='__main__':unittest.main()
