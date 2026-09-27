# Variance-component optimization for the working matched model

`scripts/fit_matched_mixed_model.py` adds optimization to the previously checked
nested background/family plus phylogenetic likelihood evaluator. It profiles
overall scale with residual component fixed to one, estimating background,
family-component and species variance ratios. It does not yet fit project data.

The implementation evaluates all eight faces defined by setting any subset of
the three ratios to zero. The all-zero face is exact; each of the seven other
faces receives three L-BFGS-B starts at ratios 0.05, 1 and 20. All 22 candidate
records are returned, including failed optimizer attempts. A failed attempt
with a better objective is retained as the best candidate and forces review,
rather than being hidden behind a successful but worse fit.

The numerical parameterization is log(1 + ratio), which admits exact zero
components. The explicit default maximum ratio is 10,000. Contact with this
bound forces review; it is not interpreted as an estimated finite optimum or
silently resolved by increasing the bound. Optimization records convergence
messages, evaluation counts, boundary flags and a separately evaluated finite
difference gradient. A lower-bound derivative must not indicate improvement
into the feasible region. All three full-face starts must agree to 1e-5 in
negative profiled restricted likelihood, and maximum absolute projected
gradient must be at most 1e-3. Passing these checks does not prove the global
optimum or parameter identifiability.

`scripts/check_matched_mixed_optimizer.py` generates three numerical fixtures:
residual-only, mixed components and strong background dependence. Each has
96 observations, nested control/family identifiers and a four-column species
factor. All 66 candidate likelihoods are independently recomputed using full
dense covariance solves and determinants. Dense generalized least-squares
coefficients and profiled scale also agree. A separate Powell optimizer with
three starts per case reaches a matching best objective; the largest absolute
difference is 2.05e-8. All three retained candidates pass the internal numerical
optimization checks. This is not a test that estimates exactly recover their
generating parameters from a finite random sample.

Complete candidate and reference outputs are in
`results/model_validation/matched-mixed-optimizer-20260927-v1`.
The checksum-bound receipt is
`metadata/matched_mixed_optimizer_checks_20260927.json`.

These computational fixtures are not a biological pilot. Complete project-data
fits, runtime planning for the full sensitivity grid, optimizer behavior on
real designs, residual diagnostics, nonlinear sequence effects, variance
identifiability and inferential calibration remain necessary. Conditional
coefficient covariance from a candidate fit must not be promoted to a final
confidence interval or significance test.
