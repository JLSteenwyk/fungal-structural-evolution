#!/usr/bin/env python3
import unittest
import numpy as np
from ancestral_state_patterns import group_patterns,reconstruct_patterns
from ancestral_categorical_diagnostics import diagnose_states


class PatternTests(unittest.TestCase):
    def test_full_roundtrip_and_duplicate_mapping(self):
        rng=np.random.default_rng(20260928)
        x=rng.integers(0,3,size=(4,30,2,6),dtype=np.uint8)
        x[:,:,1,5]=x[:,:,0,0]
        patterns,mapping=group_patterns(x)
        self.assertEqual(len(patterns),11)
        self.assertEqual(mapping[1,5],mapping[0,0])
        np.testing.assert_array_equal(reconstruct_patterns(patterns,mapping),x)

    def test_equal_frequencies_different_temporal_traces(self):
        x=np.zeros((4,40,2),dtype=np.uint8)
        x[:,20:,0]=1;x[:,::2,1]=1
        patterns,mapping=group_patterns(x)
        self.assertEqual(len(patterns),2)
        self.assertNotEqual(mapping[0],mapping[1])

    def test_chain_order_and_state_labels_preserved(self):
        x=np.zeros((4,30,3),dtype=np.uint8)
        x[0,:,0]=1;x[1,:,1]=1;x[0,:,2]=2
        self.assertEqual(len(group_patterns(x)[0]),3)

    def test_all_coordinate_diagnostics_identical(self):
        rng=np.random.default_rng(42)
        x=rng.integers(0,3,size=(4,100,8),dtype=np.uint8);x[:,:,4:]=x[:,:,:4]
        patterns,mapping=group_patterns(x)
        reports=[diagnose_states(p,'AC-') for p in patterns]
        for i in range(8):self.assertEqual(reports[mapping[i]],diagnose_states(x[:,:,i],'AC-'))

    def test_invalid(self):
        for x in [np.zeros((4,30),dtype=np.uint8),np.zeros((4,0,3),dtype=np.uint8),np.zeros((4,30,3),dtype=float)]:
            with self.assertRaises(ValueError):group_patterns(x)
        with self.assertRaises(ValueError):reconstruct_patterns(np.zeros((1,4,20),dtype=np.uint8),np.array([1]))


if __name__=='__main__':unittest.main()
