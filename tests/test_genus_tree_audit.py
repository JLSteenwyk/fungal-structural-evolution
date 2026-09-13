import io
from pathlib import Path
import sys
import unittest
from Bio import Phylo
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_genus_codon_trees import split_map


class GenusTreeAuditTests(unittest.TestCase):
    def parse(self, text):
        return Phylo.read(io.StringIO(text), 'newick')

    def test_display_root_does_not_change_split_identity(self):
        taxa = set('ABCD')
        first = split_map(self.parse('(A:1,B:1,(C:1,D:1)80/90:2);'), taxa)
        second = split_map(self.parse('(C:1,D:1,(A:1,B:1)80/90:2);'), taxa)
        self.assertEqual(set(first), {('A', 'B')})
        self.assertEqual(set(first), set(second))

    def test_malformed_or_unresolved_trees_fail(self):
        for text in ['(A:1,A:1,(C:1,D:1):2);', '(A:-1,B:1,(C:1,D:1):2);',
                     '(A:1,B:1,C:1,D:1);', '(A:1,B:1,(C:1,E:1):2);',
                     '(A,B:1,(C:1,D:1):2);']:
            with self.subTest(text=text), self.assertRaises(ValueError):
                split_map(self.parse(text), set('ABCD'))


if __name__ == '__main__':
    unittest.main()
