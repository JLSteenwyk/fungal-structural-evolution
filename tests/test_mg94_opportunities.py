import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from normalize_genus_mg94_branches import opportunities


class OpportunityTests(unittest.TestCase):
    def test_codons_preceding_stops_retain_opportunities(self):
        for code in [1,12]:
            codons,counts=opportunities(code)
            for codon,expected in [('GTT',(1,2)),('TAC',(1/3,2)),('TCT',(1,2))]:
                with self.subTest(code=code,codon=codon):
                    np.testing.assert_allclose(counts[codons.index(codon)],expected,atol=1e-12)

    def test_stop_neighbors_are_excluded_without_renormalizing_alternatives(self):
        codons,counts=opportunities(1)
        self.assertEqual(len(codons),61)
        self.assertNotIn('TAA',codons)
        np.testing.assert_allclose(counts[codons.index('TGG')],(0,7/3),atol=1e-12)
        np.testing.assert_allclose(counts[codons.index('ATG')],(0,3),atol=1e-12)

    def test_code12_ctg_reassignment_changes_opportunities(self):
        codons,counts=opportunities(1)
        np.testing.assert_allclose(counts[codons.index('CTG')],(4/3,5/3),atol=1e-12)
        codons,counts=opportunities(12)
        np.testing.assert_allclose(counts[codons.index('CTG')],(0,3),atol=1e-12)


if __name__=='__main__':unittest.main()
