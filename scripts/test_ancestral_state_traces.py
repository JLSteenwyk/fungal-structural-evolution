#!/usr/bin/env python3
import unittest
import numpy as np
from prepare_ancestral_state_traces import project_states,ALPHABET


class StateTraceTests(unittest.TestCase):
    def test_shared_anchors_gaps_and_unanchored_residue(self):
        observed={'a':'AA','b':'AA'}
        x,u=project_states({'a':'A-A','b':'-AA','n':'C-D'},observed,{'ancestor':'n'})
        self.assertEqual(''.join(ALPHABET[i] for i in x[0]),'CD-D')
        self.assertEqual(u.tolist(),[0])
        y,v=project_states({'a':'A--A','b':'--AA','other_label':'CE-D'},observed,{'ancestor':'other_label'})
        np.testing.assert_array_equal(x,y)
        self.assertEqual(v.tolist(),[1])

    def test_input_x_retains_position(self):
        observed={'a':'AX'}
        x,_=project_states({'a':'AA','n':'CD'},observed,{'ancestor':'n'})
        y,_=project_states({'a':'AG','n':'CD'},observed,{'ancestor':'n'})
        np.testing.assert_array_equal(x,y)

    def test_known_residue_change_rejected(self):
        with self.assertRaises(ValueError):
            project_states({'a':'G','n':'C'},{'a':'A'},{'ancestor':'n'})


if __name__=='__main__':unittest.main()
