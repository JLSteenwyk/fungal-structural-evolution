# Covariance dependency review

The original uniform covariance producer finished all 4,340 cohorts,
130,200 designs, 1,302,000 design/loading/tree audits and 6,220,800 setting
links. It reports every basis as `covariance_basis_requires_review`.
The independent full reader and provenance closure remain pending; a
completed producer is not an accepted model. Original jobs, numerical
tolerances, output records and failed historical stages remain unchanged.

A bounded diagnostic checked **all 300 audits in the first cohort**:
30 designs, both endpoint loading modes and five trees. The raw and
fixed-effect-projected Grams are consistent, within their recorded error
envelopes, with these kernel relations:

\[
K_{pair}=I+K_{background},\qquad
K_{gene}=\tfrac12K_{pair},\qquad
K_{model}=K_{gene}.
\]

The largest Gram-vector residual is 6.19e-11. Restricting to residual,
background, family intercept and species gives a minimum normalized Gram
singular value above 0.246 in this cohort. These are diagnostics for one
cohort, not full-grid exact operator identities or qualification of a
four-parameter model.
[Bounded diagnostic](../metadata/full_covariance_first_cohort_dependency_diagnostic_20261003_v1.json).

If the identities are independently proved for every cohort, the uniform
model can combine nonnegative variance terms:

\[
a=\theta_{pair}+\tfrac12(\theta_{gene}+\theta_{model}),\qquad
\eta_{residual}=\theta_{residual}+a,\qquad
\eta_{background}=\theta_{background}+a.
\]

Together with the already folded target/residual and family/intercept terms,
this would preserve the complete nonnegative covariance cone: every original
variance vector maps into nonnegative composite variances, and every
nonnegative composite vector is representable by setting the removed
variances to zero. It would not identify separate gene, model or model-pair
variance estimates. This is a proposed reparameterization conditional on
full identity proofs, not a model already selected or fitted.

Next steps are exact sparse-operator/cohort proofs across all 4,340 cohorts,
independent full readback, and full raw/REML basis qualification with retained
review dispositions. Additional dependencies can remain in individual
cohorts. Nonuniform residual weighting needs its own derivation; the uniform
fold must not be adopted automatically. Fitting, full optimization readback,
calibration, phylogenetic/model acceptance and biological interpretation
remain gated. This review does not complete any of the eight scientific aims.
