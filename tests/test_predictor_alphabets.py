import copy
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_predictor_alphabets import compare_states


class PredictorAlphabetTests(unittest.TestCase):
    def fixture(self):
        valid=np.ones(100,dtype=bool);valid[[0,-1]]=False
        return {'sequence':'A'*100,'states_array':np.array(list('A'*100)),
                'valid':valid,'feature_min_plddt':np.full(100,90.),
                'feature_max_pae':np.full(100,2.),'partner_residue_1based':np.full(100,50)}

    def test_identity_and_invalid_terminal_symbols(self):
        a=self.fixture();b=copy.deepcopy(a);b['states_array'][[0,-1]]='D'
        r,mask=compare_states(a,b,70,10)
        self.assertEqual(r['retained_residues'],98)
        self.assertEqual(r['state_mismatches'],0)
        self.assertEqual(r['status'],'compared')

    def test_joint_confidence_and_partner_strata(self):
        a=self.fixture();b=copy.deepcopy(a)
        b['states_array'][1:4]='C';b['partner_residue_1based'][1]=60
        b['feature_max_pae'][3]=12
        r,_=compare_states(a,b,70,10)
        self.assertEqual(r['state_mismatches'],2)
        self.assertEqual(r['same_partner_state_mismatches'],1)
        self.assertEqual(r['changed_partner_state_mismatches'],1)
        self.assertEqual(r['retained_residues'],97)

    def test_coverage_and_sequence_identity(self):
        a=self.fixture();b=copy.deepcopy(a);b['feature_min_plddt'][1:60]=20
        r,_=compare_states(a,b,70,10)
        self.assertEqual(r['status'],'insufficient_coverage')
        b['sequence']='C'*100
        with self.assertRaises(ValueError):compare_states(a,b,70,10)


if __name__=='__main__':unittest.main()
