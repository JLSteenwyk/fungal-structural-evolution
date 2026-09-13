import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_paired_marker_fits import support_options, verify_bootstrap_tips


class PairedSupportTests(unittest.TestCase):
    def test_only_searched_sequence_tree_gets_support(self):
        options = support_options('aa', 1000, 1000)
        self.assertIn('--alrt', options)
        self.assertIn('--bnni', options)
        self.assertIn('--boot-trees', options)
        for label in ('3di_af', '3di_af_empirical', '3di_llm'):
            self.assertEqual(support_options(label, 1000, 1000), [])
        self.assertEqual(support_options('aa', 0, 0), [])

    def test_invalid_replication_rejected(self):
        for a, b in ((-1, 0), (0, -1), (999, 0), (0, 999)):
            with self.assertRaises(ValueError):
                support_options('aa', a, b)

    def test_all_bootstrap_identities_and_count_required(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'trees.nwk'
            path.write_text('(A:1,B:1,(C:1,D:1):1);\n(A:1,C:1,(B:1,D:1):1);\n')
            self.assertEqual(verify_bootstrap_tips(path, set('ABCD'), 2), 2)
            with self.assertRaises(ValueError):
                verify_bootstrap_tips(path, set('ABCD'), 3)
            with self.assertRaises(ValueError):
                verify_bootstrap_tips(path, set('ABCE'), 2)
            path.write_text('(A:1,B:1,(C:1,D:1,A:1):1);\n')
            with self.assertRaises(ValueError):
                verify_bootstrap_tips(path, set('ABCD'), 1)


if __name__ == '__main__':
    unittest.main()
