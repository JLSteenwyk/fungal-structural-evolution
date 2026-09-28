#!/usr/bin/env python3
import copy
from io import StringIO
import unittest
from run_fastml_optimizer_refinement import effective_settings


class ContractTests(unittest.TestCase):
    def setUp(self):
        self.item=dict(expected_effective_options=dict(_optimizationLevel='mid',_maxNumOfIterations='100',_maxNumOfIterationsModel='100',_epsilonOptimizationModel='0.000001',_isInitGainLossByEmpiricalFreq='0'),initial_parameters=dict(alpha=1.,gain=1.,loss=1.))
        options={**self.item['expected_effective_options'],'_performOptimizationsBBL':'0','_performOptimizationsROOT':'0','_isRootFreqEQstationary':'1','_performOptimizationsManyStarts':'0'}
        self.stdout='\n'.join(k+'\t('+('Str' if k=='_optimizationLevel' else 'Float')+')\t'+v for k,v in options.items())
        self.log='### optimization starting- epsilonOptParam=1e-06 epsilonOptIter= 3e-06, MaxNumIterations=100\n L= -10 gainLossRatio= 1 gain= 1 loss= 1 Alpha= 1\n'

    def test_expected_settings_and_printed_start(self):
        self.assertEqual(effective_settings(self.stdout,self.log,self.item)['model_optimizer_calls'],1)

    def test_silent_low_level_override_rejected(self):
        with self.assertRaises(AssertionError):effective_settings(self.stdout.replace('(Str)\tmid','(Str)\tlow'),self.log,self.item)
        with self.assertRaises(AssertionError):effective_settings(self.stdout,self.log.replace('MaxNumIterations=100','MaxNumIterations=1'),self.item)

    def test_overwritten_start_and_branch_optimization_rejected(self):
        with self.assertRaises(AssertionError):effective_settings(self.stdout,self.log.replace('gain= 1','gain= 0.6'),self.item)
        with self.assertRaises(AssertionError):effective_settings(self.stdout.replace('_performOptimizationsBBL\t(Float)\t0','_performOptimizationsBBL\t(Float)\t1'),self.log,self.item)


if __name__=='__main__':unittest.main()
