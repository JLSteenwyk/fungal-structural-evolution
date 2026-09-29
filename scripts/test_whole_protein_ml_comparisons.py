#!/usr/bin/env python3
"""Manually specified contrasts protect exclusion and numerical-review behavior."""
import unittest
from export_whole_protein_ml_comparisons import comparison, PASS


def fit(nll, p=2, status=PASS, records=100):
    return dict(negative_profiled_ml=nll,fixed_coefficients=p,effective_status=status,
                records=records,support_classification='zero_supported_to_numeric_tolerance')


class Contrasts(unittest.TestCase):
    def test_nested_gain_and_added_coefficient_count(self):
        r=comparison(fit(110),fit(100,3),'left_named_columns_nested_in_right')
        self.assertEqual(r['log_likelihood_gain_right_vs_left'],10)
        self.assertEqual(r['twice_log_likelihood_gain_right_vs_left'],20)
        self.assertEqual(r['nested_model_log_likelihood_gain'],10)
        self.assertEqual(r['nested_added_fixed_coefficients'],1)
        self.assertEqual(r['comparison_status'],'numerically_checked_comparison_not_inferential_acceptance')

    def test_reverse_nesting_retains_orientation(self):
        r=comparison(fit(100,3),fit(110),'right_named_columns_nested_in_left')
        self.assertEqual(r['log_likelihood_gain_right_vs_left'],-10)
        self.assertEqual(r['nested_model_log_likelihood_gain'],10)

    def test_different_observations_have_no_direct_gain(self):
        r=comparison(fit(100),fit(10,records=50),'different_observations_no_direct_comparison')
        self.assertIsNone(r['log_likelihood_gain_right_vs_left'])
        self.assertEqual(r['comparison_status'],'different_observations_no_direct_comparison')

    def test_worse_nested_model_is_not_clamped(self):
        r=comparison(fit(100),fit(110,3),'left_named_columns_nested_in_right')
        self.assertEqual(r['nested_model_log_likelihood_gain'],-10)
        self.assertEqual(r['comparison_status'],'nested_or_identical_likelihood_violation_requires_review')

    def test_unresolved_optimizer_is_not_qualified(self):
        r=comparison(fit(110,status='ml_candidate_requires_optimization_review'),fit(100,3),'left_named_columns_nested_in_right')
        self.assertEqual(r['comparison_status'],'optimization_review_required')

    def test_missing_fit_is_explicit(self):
        r=comparison(fit(None,status='fit_error_requires_review'),fit(100,3),'left_named_columns_nested_in_right')
        self.assertIsNone(r['log_likelihood_gain_right_vs_left'])
        self.assertEqual(r['comparison_status'],'missing_fit_requires_review')

    def test_nonnested_gain_has_no_added_coefficient_contrast(self):
        r=comparison(fit(110),fit(100),'same_observations_no_named_column_nesting')
        self.assertEqual(r['log_likelihood_gain_right_vs_left'],10)
        self.assertIsNone(r['nested_model_log_likelihood_gain'])
        self.assertIsNone(r['nested_added_fixed_coefficients'])

    def test_identical_design_disagreement_is_reviewed(self):
        r=comparison(fit(110),fit(100),'identical_named_design')
        self.assertEqual(r['nesting_status'],'identical_design_likelihood_disagreement')
        self.assertEqual(r['comparison_status'],'nested_or_identical_likelihood_violation_requires_review')


if __name__=='__main__': unittest.main()
