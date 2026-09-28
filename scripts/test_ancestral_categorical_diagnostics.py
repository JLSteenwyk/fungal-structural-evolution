#!/usr/bin/env python3
import unittest
import numpy as np
from ancestral_categorical_diagnostics import diagnose_states


class CategoricalTests(unittest.TestCase):
    def setUp(self):self.rng=np.random.default_rng(20260928)

    def test_independent_categorical_and_unseen_state(self):
        x=self.rng.integers(0,3,size=(4,4000));r=diagnose_states(x,'ACD-')
        self.assertEqual(r['status'],'observed_state_indicator_screens_pass_only')
        self.assertEqual(r['unobserved_states'],['-'])
        self.assertEqual(r['indicators']['-']['screen']['status'],'state_not_observed_probability_unresolved')

    def test_constant_agreement_not_pass(self):
        r=diagnose_states(np.zeros((4,75),dtype=int),'AC-')
        self.assertEqual(r['status'],'no_observed_state_variation_requires_review')
        self.assertEqual(r['indicators']['A']['screen']['status'],'state_constant_in_all_chains_requires_review')

    def test_different_stuck_modes(self):
        x=np.repeat(np.arange(4)[:,None],100,axis=1);r=diagnose_states(x,'ACD-')
        self.assertEqual(r['status'],'categorical_mixing_requires_review')
        self.assertTrue(all(v['total_variation']==1 for v in r['pairwise_empirical_total_variation']))

    def test_frequency_shift(self):
        x=self.rng.integers(0,2,size=(4,4000));x[0]=self.rng.binomial(1,.95,size=4000)
        r=diagnose_states(x,'A-');self.assertEqual(r['status'],'categorical_mixing_requires_review')
        self.assertGreater(r['indicators']['A']['screen']['rhat'],1.01)

    def test_slow_mixing(self):
        x=np.repeat(self.rng.integers(0,2,size=(4,40)),100,axis=1);r=diagnose_states(x,'AC')
        self.assertEqual(r['status'],'categorical_mixing_requires_review')
        self.assertLess(r['indicators']['A']['screen']['bulk_ess'],400)

    def test_rare_state_missing_in_chain(self):
        x=np.zeros((4,100),dtype=int);x[0,0]=1;r=diagnose_states(x,'AC')
        self.assertEqual(r['status'],'categorical_mixing_requires_review')
        self.assertEqual(r['indicators']['C']['screen']['status'],'constant_chain_requires_review')

    def test_label_reordering_invariance(self):
        x=self.rng.integers(0,3,size=(4,2000));a=diagnose_states(x,'ACD');b=diagnose_states(2-x,'DCA')
        self.assertEqual(a['status'],b['status']);self.assertEqual(a['pairwise_empirical_total_variation'],b['pairwise_empirical_total_variation'])
        self.assertEqual(a['indicators'],b['indicators'])

    def test_invalid_and_short(self):
        for x in [np.zeros((3,100),int),np.zeros((4,0),int),np.ones((4,100),float),np.full((4,100),-1),np.full((4,100),3)]:
            with self.assertRaises(ValueError):diagnose_states(x,'ACD')
        x=self.rng.integers(0,2,size=(4,10));r=diagnose_states(x,'AC')
        self.assertEqual(r['status'],'categorical_mixing_requires_review')
        self.assertTrue(all(v['screen']['status']=='insufficient_retained_draws' for v in r['indicators'].values()))


if __name__=='__main__':unittest.main()
