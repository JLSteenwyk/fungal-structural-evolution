import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_paired_marker_fits import tree_edges


class EdgeCorrespondenceTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'tree.nwk'
            p.write_text(text)
            return tree_edges(p, set('ABCD'))

    def test_display_root_does_not_change_edge_identity_or_length(self):
        a = self.parse('(A:1,B:2,(C:3,D:4):5);')
        b = self.parse('((A:1,B:2):2,(D:4,C:3):3);')
        self.assertEqual(a, b)
        self.assertEqual(a[('A', 'B')], 5)

    def test_topology_and_invalid_lengths_are_detected(self):
        a = self.parse('(A:1,B:2,(C:3,D:4):5);')
        b = self.parse('(A:1,C:2,(B:3,D:4):5);')
        self.assertNotEqual(set(a), set(b))
        with self.assertRaises(ValueError):
            self.parse('(A:-1,B:2,(C:3,D:4):5);')


if __name__ == '__main__':
    unittest.main()
