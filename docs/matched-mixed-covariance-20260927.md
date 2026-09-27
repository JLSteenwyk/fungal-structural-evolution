# Matched-comparison covariance and restricted-likelihood evaluator

`scripts/matched_mixed_covariance.py` implements a numerical evaluator for a
working Gaussian mixed model. It does not optimize variance components or fit
the project data. The profiled restricted-likelihood framework is described by
[Bates et al. (2015), Fitting Linear Mixed-Effects Models Using lme4](https://www.jstatsoft.org/article/view/v067i01).
Our implementation specializes the covariance calculation to nested background
and family-component groups plus the previously verified species factor.

For a single fixed setting, the covariance is

\[
V=eI+bZZ^T+fHH^T+pFF^T.
\]

Z identifies reused background nodes, H identifies family components, and F
contains the chosen tree's species-factor rows for the retained observations.
e must be strictly positive; b, f and p can equal zero. Values denote variances,
not standard deviations. This working model assumes Gaussian random components
and a common residual variance. It does not establish independent residuals,
capture all shared protein/model identities within families, or validate the
additive endpoint phylogenetic assumption.

Each background must belong to one family component; the constructor rejects
violations instead of silently assigning a nesting. It first solves the
background blocks analytically, then adds the nested family intercept using
its diagonal group information, and finally adds the species term with a
Woodbury update and a small Cholesky factor. Log determinants use the
corresponding determinant identities. The full record covariance is never
constructed by the production evaluator. Zero variance components and
rank-deficient or empty species factors are supported without jitter.

At fixed component ratios, `profiled_reml` computes generalized least-squares
coefficients and profiles an overall variance multiplier using residual
quadratic form divided by n minus fixed-effect rank. It returns the negative
profiled restricted log likelihood and the conditional coefficient covariance.
The normalization contains log determinant of X-transpose V-inverse X; changing
fixed-effect bases changes a constant in the criterion. Comparing restricted
likelihoods across different fixed-effect models is therefore not justified.
The scale convention must be fixed when optimizing component ratios (for
example residual component e = 1 before profiling overall scale).

Rank-deficient fixed-effect designs, nonpositive residual variance, negative or
nonfinite component values, malformed arrays and nonpositive residual quadratic
forms fail explicitly. A conditional coefficient covariance is not a calibrated
interval accounting for estimated variance components, phylogenetic uncertainty
or multiple testing.

`scripts/check_matched_mixed_covariance.py` compares the evaluator with direct
dense covariance solves and determinants in 243 cases: three residual scales,
zero and nonzero background/family/species components, and empty, full and
duplicated species-factor columns. Every case is also checked under a random
row permutation. Five invalid-input cases must fail. Maximum absolute errors
were 9.81e-13 for solves, 3.98e-13 for log determinants and 6.26e-13 for profiled
restricted likelihood. Coefficients, conditional covariance and profiled scale
also agree with dense calculations. The receipt is
`metadata/matched_mixed_covariance_checks_20260927.json`.

These are numerical fixtures, not a reduced biological pilot. Optimizer
boundary handling, competing solutions, adequacy simulations, complete data
fits, residual dependence and uncertainty calibration remain necessary. No
inferred duplication effect, p-value or final interval is produced here.
