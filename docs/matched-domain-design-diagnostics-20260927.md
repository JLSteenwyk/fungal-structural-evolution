# Design checks before fitting matched structural contrasts

The full job covers 192 boundary, mask, cohort, coverage and input-order settings,
with 432 guide/policy/scenario strata per setting (82,944 designs). It reads the
verified domain-averaged configuration measurements and joins every qualifying
selected record. All matched-record denominators must agree with the completed
record summaries. No response values determine which designs are examined.

Four candidate nuisance covariates are examined: target-minus-background exact
sequence identity, original-interval coverage difference, log aligned-length
ratio, and joint high-confidence residue-fraction difference. Measurements are
averaged equally over eligible domains within each record, as in the verified
summaries. These are alignment-derived covariates; their availability does not
establish that they remove sequence or prediction confounding.

For each design the producer exports raw marginal minima, maxima and means;
constant-column flags; the singular spectrum after centering and scaling by
population standard deviation; rank and condition number; and whether zero lies
within every marginal covariate range. Columns with absolute range at most
1e-12 are explicitly flagged as constant. After removing those columns, the
singular values of the scaled matrix divided by square root of record count
are compared to 1e-7 for numerical rank. The intercept adds one to that rank.
Constant and deficient columns remain visible; nothing is automatically fitted
or assigned a coefficient of zero. Positive reweighting preserves exact rank,
but these equal-record condition numbers do not characterize other weightings.

The zero-range diagnostic is necessary but insufficient for support of a
contrast at equal covariates. Zero inside each marginal range need not be
inside the joint convex hull. Centering for numerical rank does not redefine
the biological equal-covariate contrast. A well-conditioned matrix does not
establish independent replication, adequate uncertainty, or model adequacy.

Each design also reports distinct targets, backgrounds, families, combined
shared-identity family components and focal taxa, plus maximum target and
background reuse. Counts describe the observations remaining after that
setting's structural qualification, and are not effective sample sizes.

Implementation: `scripts/assess_matched_domain_designs.py`. Its fixtures include
duplicated covariates, constant columns, an intercept-only design, and a full
rank random design. Full readback uses independent DuckDB joins and group
moments, with covariance eigendecomposition in place of direct SVD:
`scripts/readback_matched_domain_designs.py`. It checks every exported field.
The readback waits for the exact recorded producer process and requires a
successful terminal systemd state before reading outputs.

Plan: `metadata/matched_domain_design_plan_20260927.json`.
Launches: `metadata/matched_domain_design_launch_20260927.json` and
`metadata/matched_domain_design_readback_launch_20260927.json`.
Output: `results/structural_comparisons/matched-domain-design-diagnostics-20260927-v1`.
Expected proof: `metadata/matched_domain_design_readback_20260927.json`.
Both jobs use one CPU, 12 GiB memory, no swap and no GPU. Planning allowance is
0.1–2 hours per compute stage, with at most 1 GiB output each.

At launch the full diagnostics and independent readback are pending. They do
not provide a fitted phylogenetic effect, confidence interval, selection test
or causal duplication inference. Phylogenetic covariance and shared-control
dependence must still be specified and validated in the inferential model.
