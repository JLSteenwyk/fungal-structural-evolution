# Joint covariate support for the comparative-model intercept

The comparative models interpret the intercept at zero sequence-identity,
original-coverage, log aligned-length and confidence-fraction differences.
The completed design diagnostics place zero within every marginal covariate
range. That does not establish that the four zero differences are jointly
supported by the observed records. This analysis checks the convex hull of
record covariates for every exact production input before interpreting the
adjusted contrasts.

All 28,808 unique record inputs are included, reusing the verified inventory
that represents all 82,944 settings and five tree choices. This diagnostic
depends on covariates, so the result applies to every tree fit for the same
input. The original fits and their flags are not changed. Each reconstructed
matrix and ordered node/family/species identity sequence must match its exact
inventory fingerprint.

## Numerical definition

For record covariates X, each column is divided by its maximum absolute value;
an identically zero column has scale one. There is no centering or constant
removal. The solver minimizes t subject to nonnegative record weights summing
to one and `|sum_i w_i X_ij / scale_j| <= t` for each covariate. It uses the
installed SciPy 1.15.3 HiGHS dual-simplex interface with primal/dual feasibility
tolerances of 1e-9 and a 30-second per-input limit. See the
[official solver interface](https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html).

Every successful solve retains nonzero weights and record indices, scales,
barycenter, objective, separating direction and its minimum projection across
all records. Certificates are checked directly against the input matrix:

- A valid convex combination within 1e-8 of zero is reported as numerical
  support. This includes boundary support; it does not prove dense local data.
- A separating vector with minimum projection greater than 1e-7 establishes
  that zero lies outside the observed joint convex hull at that tolerance.
- Solver failures, invalid certificates and the intervening numerical boundary
  region remain separate unresolved dispositions. They are not omitted.

The certificates do not establish model adequacy, effective sample size,
causal exchangeability or valid uncertainty. Support is descriptive geometry of
the observed matched records, not a substitute for phylogenetic controls.

## Checks and execution

Fixtures include an interior square, a boundary triangle, an all-zero point,
a constant nonzero covariate, and the marginal-only counterexample
`(1, -0.1), (-0.1, 1), (1, 1)`. The last example contains zero in each marginal
range but has a minimum joint scaled distance of 0.45. All five expected results
and ten positive-rescaling checks passed. Synthetic input timings for 148,
2,320 and 10,963 records were approximately 0.0025, 0.017 and 0.077 seconds;
these are not full-run ETAs. Evidence:
`metadata/joint_covariate_support_checks_20260927.json`.

Run `scripts/assess_full_joint_covariate_support.py --plan
metadata/full_joint_covariate_support_plan_20260927.json`. The full run is active
with one CPU, 16 GiB memory, no swap and a 1 GiB output allowance. The planning
range is 0.25–24 hours including reconstruction and verification of the complete
source tables. The launch record pins its exact process identity. Output is
`results/model_validation/full-joint-covariate-support-20260927-v1`.

The full results, serialized-certificate readback and mapping of classifications
back to all sensitivity settings remain pending. No outcome or covariate setting
was selected because its contrast was favorable. GPU prediction remains paused.
