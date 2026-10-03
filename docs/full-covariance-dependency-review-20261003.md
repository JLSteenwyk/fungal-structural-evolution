# Covariance dependency proof and numerical qualification

Exact operator proofs have closed for all **4,340 cohorts, 75,188 logical
cases and both loading modes**: 8,680 certificates covering 68,220,240
cohort-row occurrences. Producer, independent reader, source/artifact hashes
and both original execution journals passed. Fresh observation rehashed all
4,403 bindings in the generalized proof archive. This proves a covariance
parameterization; biological model fitting and acceptance remain pending.
[Generalized proof](../metadata/full_exact_covariance_folds_completed_20261003_v2.json).
[Full binding and original-process verification](../metadata/full_exact_covariance_folds_execution_checkpoint_20261003_v2_1718.json).

## Complete-cohort identities

In every cohort and mode, the exact incidence-operator Grams satisfy

\[
K_{target}=I,\qquad
K_{gene}=\tfrac12(I+K_{background}),\qquad
K_{model}=\tfrac12K_{pair}.
\]

The signed endpoint-family kernel is zero; the unsigned kernel is exactly
four times the contrast-family-intercept kernel. The additional identity
\(K_{pair}=I+K_{background}\) holds in **3,616 cohorts per mode** and fails in
**724 per mode**. Accordingly, 7,232 certificates retain four kernels and
1,448 retain five. Every failed identity retains its original kernel and
an exact integer counterexample. All original cohorts and cases remain.

The first failed cohort is
`0048acecedff58f822f00506fdc275281324bd8f298921cd58b25dd3fc2c5c63`.
Its pair identity differs at two symmetric off-diagonal entries, including
local rows 121 and 8,385 (original case rows 68,965 and 63,394). Shared-pair
correlation explains why the initial four-kernel proposal was incomplete.
Retain the original positive-semidefinite pair kernel where needed; its
difference from \(I+K_{background}\) is not a separate PSD variance kernel.

The preserved first proof tested `gene=pair/2` and `model=gene`, retaining
seven terms for the 724 counterexample cohorts. The separately versioned V2
proof establishes the broader relationships above. Both original proof runs,
their sources, plans and outputs remain unchanged.
[First proof](../metadata/full_exact_covariance_folds_completed_20261003_v1.json).
[Earlier first-cohort diagnostic](../metadata/full_covariance_first_cohort_dependency_diagnostic_20261003_v1.json).

## Nonnegative variance cone

For the original nine nonnegative variances, define

\[
\eta_r=\theta_r+\theta_{target}+\tfrac12\theta_{gene},\quad
\eta_b=\theta_b+\tfrac12\theta_{gene},\quad
\eta_p=\theta_{pair}+\tfrac12\theta_{model}.
\]

The family-intercept composite is \(\theta_f\) in signed mode and
\(\theta_f+4\theta_{endpoint\ family}\) in unsigned mode; species variance
is unchanged. These five coefficients weight the actual residual, background,
pair, family-intercept and species kernels. Where the pair identity holds,
combine \(\eta_p\) into both \(\eta_r\) and \(\eta_b\), leaving four terms.

Every certificate contains a nonnegative forward matrix and nonnegative
right inverse whose product is exactly the identity. Every original covariance
is represented, and every nonnegative composite covariance is attainable.
Separate gene/model/target variance attribution remains unidentified.
Nonuniform weighting needs separate qualification.

## Exact and independent evidence

The producer scales dyadic weights to integers and compares sparse integer
Grams without tolerance. The reader separately forms CSC column outer products
and reconstructs the maps. All 13 original operators, all 4,340 memberships
and complete original case ordering are bound to closed operator/design
archives. Relevant files and inherited archive hashes are rechecked; broader
inherited multi-million-binding files are not freshly rehashed by this stage.

Each version passed ten synthetic scenarios, 200 covariance/right-inverse
checks, twelve malformed-record/operator rejections and complete real cohort
serialization before full execution. Producer/reader caps are two CPUs,
16 GiB, no swap, one BLAS thread, 12-GiB address space, 14,400 CPU seconds and
128 MiB per file. Large certificates and journals remain outside Git history;
compact metadata binds paths and hashes. No GPU or new charges were used.
[V2 software evidence](../metadata/exact_covariance_folds_software_validation_20261003_v2.json).

## Numerical qualification remains pending

The original uniform producer generated all 1,302,000 seven-kernel audits
and 6,220,800 setting links, reporting every basis as requiring review. Its
independent latent/QR/gesvd/SQL reader and provenance closure remain pending.
Exact identities alone cannot establish raw or REML numerical identifiability
after each fixed-effect design is accounted for.

The [complete retained-kernel workflow](full-reduced-covariance-qualification-20261003.md)
is software-qualified and queued behind that original arithmetic closure.
It will select exact principal submatrices of the closed raw/REML Grams and
unchanged error envelopes, then independently verify every matrix entry,
gesvd classification and original setting link. Unresolved norms, rank
boundaries, dependencies, non-ready outcomes and failed source audits remain
explicit. This reuses the complete original calculations and preserves every
setting.

Fitting, full optimizer readback, calibration, accepted species/reconciled gene
trees, predictor/domain/PAE controls and biological interpretation remain gated.
All eight aims and the approximately 500-fungus/25-outgroup project remain
incomplete; this proof is a prerequisite for the full sequence–structure model.
