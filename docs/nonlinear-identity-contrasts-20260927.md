# Nonlinear sequence-identity contrasts

Current status (September 27): input construction, expanded design checks and
exact-observation linking are independently verified. Full joint support and
constructive follow-up are also verified, and the full [ordinary-ML comparison](matched-ordinary-likelihood-20260927.md)
has launched. The stage descriptions below retain earlier launch checkpoints;
statements about fitting being pending describe those earlier checkpoints.

The running matched models adjust linearly for target-minus-background sequence
identity. Quadratic and cubic terms will allow sensitivity analysis of that
assumption using the same measurements and all 192 comparison settings.
This stage builds candidate inputs; it does not fit or select nonlinear models.

For each configuration and power p = 1, 2, 3, the input is the equal-domain mean
of `target_identity**p - background_identity**p`, over exactly the eligible domains
used in the original summaries. Neither a power of the identity difference nor
a difference of powers of already averaged identities is generally equivalent.
Diagnostics quantify the latter discrepancy for the quadratic term.

The producer parses the complete measurements separately with DuckDB and pandas
and compares every output row. It checks configuration identities, eligible-domain
counts and linear contrasts against the audited original summaries, then reads
back every serialized Parquet file. Source receipts, audit bindings, scripts and
input hashes are pinned in `metadata/nonlinear_identity_contrasts_plan_20260927.json`.

The service runs `scripts/build_nonlinear_identity_contrasts.py --plan
metadata/nonlinear_identity_contrasts_plan_20260927.json` with the project Python.
Resources: one CPU, 8 GiB memory, no swap; planning allowance 5–60 minutes and
2 GiB output. Process identity is recorded in
`metadata/nonlinear_identity_contrasts_launch_20260927.json`.

Status: completed with successful terminal service status and all artifact hashes
verified. The 192 settings contain 9,983,040 configuration rows, all checked by
both implementations. Completion: `metadata/nonlinear_identity_contrasts_completed_20260927.json`.
Of 1,664,432 multi-domain configuration occurrences, 1,610,736 differ by more
than 1e-10 from the quadratic shortcut using already averaged identities; the
maximum absolute discrepancy is 0.1160578783. These counts repeat configurations
across sensitivity settings and are not independent biological observations.
Existing linear fits and pinned inputs remain unchanged. Before nonlinear fitting,
expanded designs need rank/conditioning checks, renewed joint support checks,
a specified model comparison and uncertainty procedure, and compute estimates.
Restricted-likelihood objectives with different fixed-effect spaces must not
be treated as directly comparable likelihood-ratio evidence. These inputs alone
establish no biological effect.

## Expanded model design checks

`assess_nonlinear_matched_designs.py` is running under the pinned plan
`metadata/nonlinear_matched_design_plan_20260927.json`. It covers both polynomial
degrees across all 82,944 settings (165,888 design rows), reports explicit constant
columns and scaled conditioning, and compares direct SVD with pivoted QR followed
by SVD. Constant, redundant-column and full-rank fixtures passed. The producer
receipt will still require independent source/output readback before acceptance.
Resources: one CPU, 12 GiB memory, no swap, 1 GiB output allowance; estimated
0.25–4 hours. Joint nonlinear support and fitting remain pending.

The independent readback is queued as
`fungal-nonlinear-matched-design-readback-20260927.service`, using one CPU and
12 GiB, with a 0.25–4 hour estimate after producer completion. It checks the exact
producer identity and terminal success before rebuilding all 165,888 rows using
SQL joins, covariance moments and eigenvalues. Full counts, reuse maxima, ranges,
means, spectra, ranks, conditions and marginal support are compared. Polynomial
and redundant-column numerical fixtures passed. Near-threshold spectral
disagreement will stop verification; it will not silently approve a design.

## Expanded input inventory running

The earlier linear-model reuse map is insufficient to establish equality of
nonlinear inputs. `inventory_nonlinear_model_inputs.py` now reconstructs every
quadratic/cubic setting, preserving the exact ordered response/covariate values
and target, background, family-component and species-pattern identities. The
polynomial degree is part of the fingerprint. Only signed zeros are canonicalized;
no tolerance-based merging is allowed. All 165,888 settings remain represented.

An independent SQL reconstruction, `readback_nonlinear_model_inputs.py`, waits
for producer terminal success and will verify all fingerprints, source labels,
record counts and representative recipes. Plans/launches use the
`metadata/nonlinear_model_input_*_20260927.json` prefix. Each stage uses one CPU
and 12 GiB, no swap, with 0.1–2 hours estimated; inventory output allowance is
2 GiB. Candidate tree-fit counts will estimate workload only: nonlinear model
fitting, joint support checks and inferential comparisons are not launched.

Ordinary ML objective and analytic-score evaluators now pass independent dense
checks; see [likelihood methods](matched-ordinary-likelihood-20260927.md). All
candidate fixed-effect spaces, including the linear reference, will require ML
fitting before likelihood comparisons. No nonlinear fit has been launched.

The design producer has now completed all 165,888 rows with no rank deficiencies
after constant removal and no zero reference outside any marginal range. These
are producer results pending the active full independent readback, and do not
establish joint nonlinear support.

## Nonlinear joint-reference support queued

The inventory producer finished successfully: 165,888 settings map to 57,616
distinct quadratic/cubic inputs, implying 288,080 candidate tree fits across the
five trees. The full independent inventory readback remains active; these are
not yet accepted production recipes. The ordinary-ML linear reference would
require additional fits for fair fixed-effect comparisons.

`assess_nonlinear_joint_support.py` waits for the exact inventory-auditor process
and successful terminal service before checking every one of the 57,616 expanded
inputs. It reconstructs original ordered identities and numeric fingerprints,
then applies the existing zero-reference convex-hull method to all covariates,
retaining constants. It records sparse support weights or separating directions
where resolved and preserves every unresolved case.

`readback_nonlinear_joint_support.py` then independently checks every saved
certificate against its reconstructed matrix and maps classifications back to
all 165,888 settings. Producer/readback use one CPU and 12 GiB each, no swap,
2 GiB output allowances; planning estimates are 1–12 hours for support and
0.25–4 hours for readback after prerequisites. Plans and exact process identities
are under `metadata/nonlinear_joint_support_*_20260927.json`.

Neither marginal range coverage nor full column rank proves joint support. Even
a valid support certificate establishes numerical convex-hull inclusion only,
not interior overlap, dense local observations, causal exchangeability or valid
confidence intervals. Unsupported cases must remain visible in model reporting.

## Inventory independently verified and linked to linear references

Both inventory services terminated successfully. Full independent SQL readback
verified all 165,888 setting rows and 57,616 distinct quadratic/cubic inputs.
`link_polynomial_model_inputs.py` then joined every setting to its linear reference,
checking shared source bindings, exact record counts, ordered observation
identities and nested base-column names. All 82,944 settings are represented by
28,808 unique linear/quadratic/cubic input triplets; no linear input splits into
multiple expanded triplets in this dataset. This was checked, not assumed from
the earlier inventory. Both exported tables passed serialized readback.

Artifacts: `results/model_validation/polynomial-model-input-links-20260927-v1/`
contains the full setting map and unique comparison sets. Completion and source
hashes: `metadata/polynomial_model_input_links_completed_20260927.json`. The
nonlinear support producer is now able to proceed through its verified dependency.
This does not authorize comparing existing REML objectives with future ML fits.

## Expanded design verification completed

Both producer and independent readback terminated successfully. All 165,888
design rows were verified, with no rank deficiencies after constant removal.
`record_nonlinear_design_completion.py` records source/terminal bindings and
full-table numerical summaries in
`metadata/nonlinear_matched_design_completed_20260927.json`.

| Identity degree | Settings | Median scaled condition | 95th percentile | Maximum | Minimum scaled singular value |
|---|---:|---:|---:|---:|---:|
| Quadratic | 82,944 | 11.56 | 13.80 | 15.60 | 0.09081 |
| Cubic | 82,944 | 110.74 | 176.54 | 218.95 | 0.007922 |

The cubic-to-quadratic condition ratio has median 9.64 across paired settings.
Adding the cubic term preserves full rank but increases predictor dependence.
These are centered, population-SD-scaled fixed-design diagnostics before
covariance weighting, not fitted coefficient uncertainty or evidence that a
polynomial degree is biologically preferred. Repeated sensitivity settings are
not independent samples. Cubic terms remain in the planned sensitivity analysis;
model comparisons and uncertainty must assess their contribution. Joint support
verification is still running and is not implied by marginal range coverage.

## Full joint support and constructive follow-up verified

All four stages terminated successfully: full support production, independent
matrix/certificate readback, failure census and constructive follow-up. Original
classifications are preserved for every input and all 165,888 settings.

| Degree | Original valid support | Original unresolved | Separate nonnegative certificates verified |
|---|---:|---:|---:|
| Quadratic | 28,744 | 64 | 64 |
| Cubic | 28,760 | 48 | 48 |

All 112 original failures involved negative candidate weights, with maximum total
negative mass 9.881e-10. The follow-up removed negative candidate weights,
renormalized the positive weights and recomputed the barycenter against the
original fingerprinted matrix. Every projected certificate passed explicit
nonnegativity, normalization, scalar/vector barycenter and serialized-output
checks; maximum scaled zero-reference residual was 5.995e-10, below 1e-8.
These 112 inputs map to 180 setting/degree rows. Original classifications remain
visible alongside the separate constructive-certificate field.

Thus all 57,616 expanded inputs have numerical zero-reference support from either
the original or a separately constructed nonnegative certificate. This does not
prove an LP optimum, interior overlap, dense local sampling, causal comparability,
adequate covariance or valid confidence intervals. Ordinary-ML fitting, review,
model adequacy and inferential calibration remain open.

`record_nonlinear_support_completion.py` verifies terminal states, every declared
artifact, receipt bindings, complete counts and exact failure/projected input
identity. Completion evidence: `metadata/nonlinear_joint_support_completed_20260927.json`.
# Interim optimization-flag census (September 28)

A fixed prefix of 75,205 completed ordinary-ML dispositions contains 601
optimization-review flags. Every flagged file and payload hash was verified;
the source manifest prefix was saved and checked for stability while the live
run continued. This is an interim census, not a full-grid audit.

The 601 flags comprise 587 unsuccessful selected optimizer results, six failed
projected-gradient checks, and eight disagreements among full-face starts.
There were no upper-bound-contact flags in this prefix. Flagged polynomial
degrees were 163 linear, 137 quadratic and 301 cubic. All 601 had a successful
alternative candidate within 1e-7 of the lowest objective, but that observation
does not resolve convergence or start-agreement concerns. Original statuses
and numerical tolerances remain unchanged. These results motivate targeted
optimizer follow-up after the full output audit; they are not evidence for a
biological sequence–structure relationship.

Script: `scripts/census_polynomial_optimization_flags.py`; frozen prefix,
per-fit reasons and hashes: `results/model_validation/polynomial-optimization-flag-census-20260928-v1`.
The reason counts were also independently recounted from each flagged payload;
evidence: `metadata/polynomial_optimization_flag_census_completed_20260928.json`.

A separate ordinary-ML refinement routine is now available as
`scripts/refine_matched_ml_analytic.py`. It retains the original unoptimized
reference, the all-zero face, the existing three starts on every nonempty
variance-component face, and one additional full-face start from the original
parameters (24 candidates total). It keeps the same variance bounds and review
criteria, increases the iteration/line-search allowance, and tightens optimizer
stopping tolerances. Every candidate is checked against the direct ML evaluator.
An exactly tied successful candidate takes precedence, but a lower unsuccessful
candidate remains selected and flagged; no approximate tie tolerance is used.

`scripts/check_matched_ml_refinement.py` passed six synthetic fits across three
polynomial degrees and zero/nonzero species factors. All 144 candidate
likelihoods matched a separate dense covariance calculation, with maximum
absolute error 2.843e-14; four invalid starting vectors were rejected. Evidence:
`metadata/matched_ml_refinement_checks_20260928.json`. No production refinement
has been launched. Full-grid audit, exact input reconstruction, resource
estimation and independent refinement-output checks are still needed; the
synthetic tests do not resolve any of the 601 observed production flags.

Exact input reconstruction for that frozen snapshot is now running separately:
`scripts/prepare_flagged_polynomial_inputs.py`. The 601 tree fits correspond
to 449 unique inputs across 25 partitions and 1,081,570 record occurrences.
Each reconstruction must match the original audited recipe's numerical-byte
and ordered-identity hashes. Exported NPZ arrays contain the response/design
matrix, background/family codes, species-factor rows and row identities, and
are read back in full. Resource limits are one CPU, 16 GiB RAM, no swap and a
2 GiB output allowance; planned runtime is 0.02–2 hours. Plan and exact launch
records use `metadata/flagged_polynomial_input_preparation_*_20260928.json`.
This prepares a frozen subset for future review, not a replacement full-grid
audit. It launches no model fits, changes no source outputs, and still requires
independent input validation before refinement.

Reconstruction completed successfully at 23:43:44 EDT: all 449 original input
hashes matched, covering 1,081,570 record occurrences and 601 flagged tree fits.
All exported artifact hashes were verified after successful process termination;
evidence: `metadata/flagged_polynomial_inputs_producer_completed_20260928.json`.
The separate checker `scripts/readback_flagged_polynomial_inputs.py` is running.
It reconstructs source membership and ordered identities, independently codes
background/family labels, verifies species-factor row mappings, compares every
numerical array with its audited recipe hash, and retains all original flagged
fit links. It is limited to one CPU and 16 GiB RAM with no swap, with 0.02–2 hours
planned. Plan and launch records use
`metadata/flagged_polynomial_input_readback_*_20260928.json`. Production fits
remain unchanged and no review flag has been resolved.

The first input checker stopped before mapping validation because it opened
the gzip-compressed pair table as plain text (`UnicodeDecodeError`). A separate
v2 checker uses the gzip text reader; all 52,675 pair rows were parsed in its
prelaunch check. The failed script and launch remain preserved. The replacement
is `scripts/readback_flagged_polynomial_inputs_v2.py`, with plan and launch
suffixes `20260928_v2.json`. No data, hashes or validation tolerances changed.

The corrected audit finished successfully at 23:46:08 EDT, verifying all 449
inputs, 1,081,570 record occurrences and 601 flagged fit links. Proof:
`metadata/flagged_polynomial_inputs_completed_readback_20260928_v2.json`.

With the frozen inputs independently validated, separate numerical follow-up
can proceed before the rest of the production grid finishes. It remains
ineligible for integration until full production audit and independent
refinement-output validation. `scripts/refine_flagged_polynomial_snapshot.py`
now runs the tested refinement on all 601 frozen flagged fits, first replaying
each original parameter vector and checking its likelihood, coefficient and
scale against the saved original. All original outputs remain unchanged.

The 601 original fits recorded 3,504.71 solver seconds (0.974 serial hours).
Refinement has a deliberately broad 1–16 hour planning allowance on one CPU,
16 GiB RAM, no swap and 1 GiB output. This is not a measured refinement ETA.
Each of the 24 candidate parameter vectors per fit receives a direct-likelihood
check; unsuccessful outcomes remain explicit. Output:
`results/model_validation/flagged-polynomial-refinements-20260928-v1`.
Plan and exact launch records use
`metadata/flagged_polynomial_snapshot_refinement_*_20260928.json`.

A separate output replay is queued behind successful refinement termination:
`scripts/audit_flagged_polynomial_refinements.py`. It replays all 14,424 candidate
likelihoods across 601 fits, checks the complete face/start grid and preserved
original parameters, recomputes gradients and review classifications, and
checks fitted coefficients, scales, covariance matrices and raw-unit
conversions. It uses the already validated direct-likelihood and gradient
libraries, so this is a separate execution and output audit, not an independent
mathematical implementation. Resource limits are one CPU, 16 GiB RAM, no swap,
and 0.1–8 planned hours after the producer. Plan and launch records use
`metadata/flagged_polynomial_refinement_readback_*_20260928.json`. No refined
production status is accepted solely from the live progress log.

## Direct gradient diagnostic for an unresolved refinement

For input `2212faa1575538e761e1aa1e1a9ff90ad892c05f1b0b75b40c623efba7f99264`
under `mafft_guide` (8,201 records), refinement still failed the projected
gradient threshold. The optimizer reported relative-objective convergence, but
the species derivative was 0.00117371047, above the fixed 0.001 criterion.
`scripts/diagnose_refinement_gradient.py` evaluated the direct likelihood at
four central-difference steps per component, preserving every estimate. At
step 1e-6 the species derivative was 0.00117324817, supporting the analytic
flag. Larger steps show truncation effects and the smallest step shows
cancellation sensitivity; none is used to relax acceptance or replace the
fit. All 25 likelihood evaluations completed successfully. Evidence:
`metadata/refinement_gradient_diagnostic_completed_20260928.json`; full table:
`results/model_validation/refinement-gradient-diagnostic-20260928-v1`.
Further numerical refinement remains necessary for this candidate.

The separate curvature diagnostic `scripts/diagnose_refinement_curvature.py`
then evaluated two finite-difference Hessians (steps 1e-5 and 1e-6), both
positive definite, and eight bounded Newton proposals. Both full steps changed
log1p variance ratios by approximately (4.265e-8, 8.504e-9, -3.620e-8).
The direct objective decreased by 2.274e-11 and maximum analytic gradient fell
to 2.836e-9 and 3.616e-9, respectively. This is evidence consistent with
premature relative-objective stopping near a stationary point; the tiny
objective improvement has no biological interpretation. Proposals remain
diagnostic: independent finite-difference confirmation and a reusable audited
acceptance procedure are needed before replacing any fitted result. All eight
trials and source hashes are preserved in
`results/model_validation/refinement-curvature-diagnostic-20260928-v1`;
completion evidence is `metadata/refinement_curvature_diagnostic_completed_20260928.json`.

Both full-step proposals subsequently passed direct finite-difference checks
at steps 1e-6 and 1e-7 in all three parameter directions (27 direct objective
evaluations including the original and proposed points). The largest absolute
estimated derivative was 1.137e-5, below the unchanged 1e-3 criterion. Both
replayed objectives were 2986.8187322513263, slightly below the original.
Checker: `scripts/check_refinement_curvature_proposals.py`; proof:
`metadata/refinement_curvature_proposals_checked_20260928.json`.
This verifies the two fixed proposals locally. Incorporation into a reusable,
audited refinement procedure and full-grid integration remain outstanding;
the original stored review flag has not been changed.

The diagnostic is now implemented in `scripts/interior_gradient_polish.py`.
It explicitly declines points near fixed bounds, rejects nonpositive curvature,
retains both Hessian-step proposals, and requires each to avoid worsening the
direct objective and pass analytic and both direct finite-difference checks
under the existing 1e-3 criterion. It does not select or replace a production
fit. Mathematical fixtures checked a known quadratic minimum, a stationary
no-op, boundary handling, invalid parameters, negative curvature, and an
intentionally inconsistent gradient rejected by direct differences. Evidence:
`metadata/interior_gradient_polish_fixture_checks_20260929.json`.

`scripts/check_interior_polish_case.py` reproduced both previously checked
real-case proposals through this guarded implementation. Its complete inputs,
output and source hashes are recorded in
`metadata/interior_gradient_polish_case_checked_20260929.json`. Application to
other flags and any integration remain gated on the full refinement audit;
boundary cases and disagreements among starts are not resolved by this local
routine.

### Audited gradient-only proposals queued (September 29)

`scripts/prepare_audited_gradient_proposals.py` is queued behind successful
termination of the exact recorded refinement-output audit. It requires the
complete 601-fit audit proof and matching source receipt before doing numerical
work. The pinned plan and live process identity are recorded in
`metadata/audited_gradient_proposals_plan_20260929.json` and
`metadata/audited_gradient_proposals_launch_20260929.json`.

Every audited fit receives a disposition. Fits already passing are retained;
other review flags are preserved. Only fits whose sole failing criterion is
the projected gradient receive proposals from the validated interior correction
routine. Both Hessian-step proposals must satisfy the unchanged gradient
threshold and direct finite-difference checks. Source and output hashes are
checked before and after execution. No proposal is selected or substituted for
an existing result. A separate proposal-output readback remains required before
integration, as does completion and audit of the original full model grid.

Resources: one CPU, 16 GiB memory, no swap, no GPU or paid resources; planning
allowance 0.1–8 hours and 0.1 GiB output. The job is currently waiting on its
upstream audit; a launch is not evidence of completed proposals.
