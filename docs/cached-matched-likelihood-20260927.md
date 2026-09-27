# Cached likelihood computation

`scripts/cached_matched_likelihood.py` implements the same working covariance
and profiled restricted likelihood as the previously checked direct evaluator,
with residual ratio fixed at one. It caches cross-products of the species
factor, fixed-effect design and response; background-group sums; cross-products
of those sums grouped by control reuse count; and a sparse background-to-family
incidence matrix. Variance-ratio evaluations combine these sufficient statistics
instead of repeatedly aggregating every original observation.

For M = [F, X, y], the background-adjusted cross-product starts at M-transpose M
and subtracts each group outer product weighted by b/(1 + b times group size).
Weighted group sums then give the nested family correction. A small Cholesky
solve applies the species correction, yielding the generalized least-squares
information, response cross-product and profiled residual quadratic form.
The determinant calculation is unchanged. Zero variance components remain
allowed; no covariance jitter or model approximation is introduced.

Subtracting explained from total response quadratic form can lose precision
near an exact fit. When the residual quadratic form is at most 1e-9 times the
sum of their magnitudes, the implementation uses the original direct residual
calculation instead. Nonfinite intermediates fail explicitly. This safeguard
does not replace conditioning diagnostics or guarantee accuracy on every input.

`scripts/check_cached_matched_likelihood.py` checks 81 cases combining empty,
ordinary and duplicated species factors with all zero/moderate/10,000-ratio
component combinations. Likelihood, coefficients, scale, residual quadratic form
and conditional coefficient covariance agree with the verified direct evaluator;
maximum observed likelihood discrepancy is below 9e-11. A large fitted-response
fixture triggers the direct fallback and reproduces its results exactly.

Five paired timing/comparison repetitions also cover numerical problems with
148, 2,320 and 10,963 records and 242 species-factor columns. Setup is measured
separately. Initial measurements show approximately 6-fold and 12-fold lower
evaluation time at the two larger sizes. Exact current timings and all source
hashes are in `metadata/cached_matched_likelihood_checks_20260927.json`; timings
vary with concurrent workloads and do not imply the same optimizer speedup.

These are computational fixtures, not a biological pilot. The existing
optimizer remains unchanged and checksum-bound to its completed checks and
running inventory. A separate cached optimizer still needs full objective,
boundary and solution-equivalence checks before production adoption. Actual
fitting, covariance adequacy, calibration and final scientific inference remain
incomplete.
