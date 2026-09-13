import io
import sys
import unittest
from pathlib import Path
from Bio import Phylo
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_marker_tree_support import supported_splits


class TreeSupport(unittest.TestCase):
    def parse(self, text):
        return Phylo.read(io.StringIO(text), 'newick')

    def test_split_identity_is_independent_of_unrooted_representation(self):
        a = supported_splits(self.parse('((A:1,B:1)90:0.1,C:1,D:1);'))
        b = supported_splits(self.parse('(A:1,B:1,(C:1,D:1)90:0.1);'))
        self.assertEqual(a, b)
        c = supported_splits(self.parse('(A:1,B:1,(C:1,E:1)90:0.1);'))
        self.assertNotEqual(a[0]['split_sha256'], c[0]['split_sha256'])

    def test_missing_support_is_not_zero_and_invalid_values_fail(self):
        missing = supported_splits(self.parse('((A:1,B:1):0,C:1,D:1);'))[0]
        zero = supported_splits(self.parse('((A:1,B:1)0:0,C:1,D:1);'))[0]
        self.assertEqual((missing['support_status'], missing['sh_alrt_percent']), ('not_reported', ''))
        self.assertEqual((zero['support_status'], zero['sh_alrt_percent']), ('reported', 0))
        for text in ['((A:1,B:1)101:0,C:1,D:1);', '((A:1,B:1)90:-1,C:1,D:1);']:
            with self.assertRaises(ValueError):
                supported_splits(self.parse(text))


if __name__ == '__main__':
    unittest.main()
