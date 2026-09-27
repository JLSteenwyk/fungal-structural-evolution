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

Status: launched; completion requires terminal service success and the full
receipt in `results/model_validation/nonlinear-identity-contrasts-20260927-v1`.
Existing linear fits and pinned inputs remain unchanged. Before nonlinear fitting,
expanded designs need rank/conditioning checks, renewed joint support checks,
a specified model comparison and uncertainty procedure, and compute estimates.
Restricted-likelihood objectives with different fixed-effect spaces must not
be treated as directly comparable likelihood-ratio evidence. These inputs alone
establish no biological effect.
