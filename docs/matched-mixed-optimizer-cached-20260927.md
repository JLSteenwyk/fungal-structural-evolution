# Cached optimizer equivalence

`scripts/fit_matched_mixed_model_cached.py` uses the checked cached likelihood
while preserving the original optimizer's eight zero-component faces, three
starts per nonempty face, parameterization, bounds, convergence settings and
review rules. The original optimizer remains unchanged. The new version also
records elapsed time and counts of cached versus direct-fallback evaluations.

`scripts/check_matched_mixed_optimizer_cached.py` repeats the original three
synthetic cases with the same seed and data generation. All 66 candidate
likelihoods are checked by independent dense covariance calculations. The
best solution in each case is compared with three independent dense Powell
starts and with the saved, checksum-verified original optimizer result.
Best objectives agree within 1e-6 absolute tolerance, coefficients within
1e-4 relative/1e-5 absolute tolerance, and zero-boundary flags and diagnostic
statuses agree exactly. The complete face/start enumeration is also identical.

All three cases passed. Differences from the best dense Powell objective are
at most 2.05e-8. Cached evaluation counts were 398, 766 and 818, with total
optimization times about 0.064, 0.134 and 0.138 seconds for these 96-record,
four-factor-column computational fixtures. These timings do not forecast
production fits with up to 242 species-factor columns. Small floating-point
differences can change evaluation counts without changing the accepted
solution; bit-identical optimizer trajectories are not claimed.

Complete outputs are in
`results/model_validation/matched-mixed-optimizer-cached-20260927-v1`.
The verification receipt is
`metadata/matched_mixed_optimizer_cached_checks_20260927.json`.
All numerical-library threads were restricted to one during validation.

This completes numerical integration of the cached evaluator. The full
record-input reuse map still requires its independent check before production
reuse. Production scheduling must retain all settings and tree alternatives,
preserve unconverged/boundary-review cases, and verify fitted objectives against
the direct evaluator. These tests do not establish biological model adequacy,
variance identifiability, calibrated uncertainty or an evolutionary effect.
