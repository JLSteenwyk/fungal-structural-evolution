# Nonlinear sequence-identity contrasts

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
