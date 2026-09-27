# Ordinary likelihood for fixed-effect comparisons

The running matched models use REML. Their objective values must not be compared
across different fixed-effect spaces to select a sequence-identity polynomial.
New, separate modules implement ordinary profiled Gaussian maximum likelihood
without changing any running script or pinned input.

For covariance `V = s² R`, where `R = I + b ZZᵀ + f GGᵀ + p FFᵀ`, the coefficients
at fixed ratios minimize `q = (y-X beta)ᵀ R⁻¹ (y-X beta)`. The ordinary ML variance
estimate is `q/n`, and the negative profiled log likelihood is

`0.5 [log|R| + n (1 + log(2 pi q/n))]`.

It has no fixed-design determinant term. For a ratio with derivative matrix `K`,
its score is `0.5 [tr(R⁻¹ K) - n uᵀ K u/q]`, with `u = R⁻¹ (y-X beta)`.
The analytic implementation differentiates cached cross-products; validation
instead constructs dense covariance matrices and uses this trace identity.
Coefficients at fixed ratios are shared with REML, so the cancellation fallback
reuses the existing direct residual calculation and explicitly changes the scale
denominator and objective. Returned fields distinguish `negative_profiled_ml`
from REML and record the variance-profile denominator separately from residual
degrees of freedom.

Implementation:

- `scripts/cached_matched_ml.py`: cached likelihood and direct residual fallback.
- `scripts/matched_ml_gradient.py`: analytic ratio and log1p-ratio derivatives.
- `scripts/check_matched_ml.py`: independent dense Gaussian tests.

Validation passed 243 cases across linear/quadratic/cubic designs, zero/moderate/
large covariance ratios and empty/duplicated species factors. All 729 component
derivatives agreed with dense scores; maximum absolute errors were 2.49e-10 for
likelihood and 7.00e-10 for derivatives. Another 243 checks verified invariance
to invertible fixed-design reparameterization, and 162 verified that adding fixed
effects cannot worsen ordinary likelihood at the same covariance. A large-mean
cancellation case used the direct residual fallback; a rank-deficient design was
rejected. Evidence: `metadata/matched_ml_checks_20260927.json`.

These are numerical evaluators, not completed biological fits. Later comparisons
must fit **all** candidate fixed-effect spaces, including the linear reference,
under ordinary ML on identical observations with the same covariance options.
The existing REML estimate cannot substitute for the ML optimum. Boundary and
multistart optimizer checks, nonlinear joint support, model adequacy, uncertainty
and multiplicity calibration remain required. No automatic chi-square calibration,
selected polynomial degree or causal duplication effect is established here.

## Analytic optimizer numerical validation

`fit_matched_ml.py` now enumerates all eight zero-component faces with three
starts on each nonempty face (22 candidates), using analytic log1p-ratio scores.
Every candidate is reevaluated with the separate cached objective. All attempts
remain recorded, including failed attempts. Acceptance requires optimizer success,
projected gradient at most 1e-3, agreement of the three full-face objectives within
1e-5, and no upper-bound contact. Bounds are never silently expanded. Analytic
cancellation fails explicitly; it does not substitute a noisy finite difference.

`check_matched_ml_optimizer.py` passed six synthetic fits across three nested
fixed designs and empty/nonempty species factors. All 132 candidate objectives
were independently checked against dense Gaussian calculations, as were final
coefficients, scale, conditional covariance and analytic scores. A balanced
random-intercept fit matched its closed-form ML solution; an exact zero-group-mean
case reached the zero-variance boundary. Forced iteration limits retained failed
attempts and required review; a deliberately restrictive variance bound also
required review. Evidence: `metadata/matched_ml_optimizer_checks_20260927.json`.

Numerical convergence does not establish identifiable variance components. The
empty-factor cases can return arbitrary species ratios on a flat likelihood;
`species_kernel_is_zero` records this explicitly and
`component_identifiability_assessed` remains false. Such a ratio is not evidence
for phylogenetic variance. Production-size timing, expanded joint support, full
input verification and statistical calibration remain before production ML fitting.

## Complete-fit size benchmark

The synthetic size benchmark finished successfully on one CPU with 138.2 MiB
peak service memory. Every case used seven fixed-effect columns, a 242-column
species factor and all 22 optimizer attempts. Complete fitting took 2.03 seconds
for 148 records, 3.66 seconds for 2,320 records, and 12.45 seconds for 10,963 records.
All three passed numerical optimization checks. Every one of the 66 candidate
objectives was separately checked using direct residual evaluation; maximum
absolute objective error was 7.28e-12. Validation took an additional 0.04, 0.47
and 2.11 seconds respectively. Source/output hashes and terminal success were
verified in `metadata/matched_ml_size_benchmark_completed_20260927.json`.

These are three synthetic timing observations, not a representative convergence
distribution or full-grid ETA. Full workload sizing still needs the completed
expanded inventory, observed record-size distribution, shared machine load and
allowance for difficult or flagged fits. The existing REML run remains unchanged.
