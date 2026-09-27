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
