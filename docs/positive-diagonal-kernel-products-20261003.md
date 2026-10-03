# Fresh positive-diagonal raw and REML kernel products

The independent numerical audit supports finite positive residual diagonals,
with target identity covariance retained separately. Software qualification
passed dense error-contrast, component-kernel and latent-overlap comparisons.
A separate reuse implementation also passed, sharing freshly computed unchanged
products across diagonals. **Full real-data weighted raw/REML qualification
and fitting remain pending.**

The prerequisite [original-cohort weight export](inverse-reuse-weight-controls-20261003.md)
has independently checked all 4,340 cohorts, 75,188 logical cases and
136,440,480 case/control occurrences. Closure includes 13,071 hashes and
two actual original journals. Working diagonals do not estimate calibrated
measurement precision or posterior ESS.

## Computation

For fixed design `X`, let `Q` span its column space orthonormally and
`H = I - Q Q'`. Covariance Frobenius products are:

```text
G_ab       = trace(K_a K_b)
G_ab(REML) = G_ab - 2 trace(Q' K_a K_b Q)
                    + trace((Q' K_a Q)(Q' K_b Q))
```

For residual kernel `D`, retain `D Q` and `Q' D Q`. Its raw norm is
`sum(D_i^2)`; uniform shortcuts `n` and `n-p` do not apply to arbitrary
diagonals. For entity factor `Z`, diagonal/entity products are
`sum_i D_i * sum_j Z_ij^2`, and entity/entity products can independently use
`||Z_a' Z_b||_F^2`. Species factors use the same latent identities.

The primary reference constructs component-local kernels. The independent
reader uses sparse latent overlaps, extended-precision design scaling and
pivoted QR with every diagonal correction. Explicit dense complete QR supplies
an error-contrast reference for software checks. Scalable implementations
avoid full global covariance/projector matrices.

Reuse computes nonresidual products and fixed-design images/cores once,
then rebuilds residual rows, projected Grams and every envelope for each `D`.
All products are fresh. Constructor ones are a disposable placeholder for
a residual row that is replaced completely; no saved uniform audit, envelope
or qualification is consumed. Source factors are copied before read-only
caching, and returned lists do not expose cached arrays. Whole-data speedup
has not been measured.

No basis is deleted on a rank result. Constant and near-uniform diagonals,
zero kernels, unresolved norms and boundaries retain review states. Any named
uniform fold requires its separate closed source certificate; nonconstant
controls require the positive-diagonal cone. Arithmetic checks select no
biological variance model.

## Qualification and evidence

The general checker covers 24 cases: signed/unsigned, five/six kernels and
six patterns—uniform one, constant two, alternating, positive ramp,
near-uniform and a `1e-4` to `1e4` diagonal. It passed 24 component/latent/dense
raw/REML comparisons, 24 permutations and 24 full-rank design transformations.
Eight constant-diagonal and four explicit-zero reviews were retained;
60 invalid diagonal/design/factor/component/name inputs were rejected.

Maximum raw dense error was 1.46e-11 and projected error 2.91e-11, below the
combined fresh roundoff envelopes. Maximum error/bound fractions were
0.000241 raw and 0.000159 projected. Tolerances remain `rtol=3e-9`,
`atol=2e-8`. Twenty-four of 48 raw/REML diagnostic outcomes are dependence
reviews, including near-uniform cases. Nothing is clipped or deleted.
Roundoff envelopes are numerical review policies, not statistical uncertainty.

Reuse passed all 24 complete separate-audit comparisons, independent latent
and dense checks, and 24 reordered calls. Eight checks cover source-factor
copies and detached returned results; 40 invalid diagonal/design cases were
rejected. Dense errors and rank/dependency reviews were unchanged.

Both actual original software waits exited zero with matching exact wrapper
PID/create/command/invocation messages and native receipts. Each used two
CPUs/16 GiB/no swap, one BLAS thread, 12-GiB native address-space,
900-CPU-second/256-MiB-per-file limits and a 1,800-second wall cap. Child CPU
times were 1.62 and 1.31 seconds; child peak RSS was 116,391,936 and
114,819,072 bytes. Different workloads are not a production speed comparison.
Manager memory accounting is distinct from child RSS.

The first general-software resource file mistakenly estimated at most 12
synthetic rows. All 24 actual cases have 40 rows; the estimate was not an
enforced input limit. A scope audit preserves the original and records 40.
Executed cases, tolerances and budgets are unchanged. The reuse resource
file records 40 rows. No launched source/proof was rewritten.

- [General validation](../metadata/positive_diagonal_kernel_software_validation_20261003_v1.json),
  [original wait/journal proof](../metadata/positive_diagonal_kernel_software_transport_20261003_v1.json)
  and [preserved size correction](../metadata/positive_diagonal_kernel_software_scope_audit_20261003_v1.json).
- [Reuse validation](../metadata/positive_diagonal_basis_context_software_validation_20261003_v1.json)
  and [original wait/journal proof](../metadata/positive_diagonal_basis_context_software_transport_20261003_v1.json).
- [Full weight completion](../metadata/full_inverse_reuse_weights_completed_20261003_v1.json)
  and [fresh whole closure/terminal verification](../metadata/full_inverse_reuse_weight_completed_verification_20261003_v1.json).

Reproduce the two `check_positive_diagonal_*` scripts under saved wrapper/
resource commands with fresh output/receipt paths and pinned source versions.
Use the existing shared-entity environment. The new modules preserve the
frozen uniform reader, failed outputs and original queues; they do not explain
the old mismatch.

Full weighted source/design linkage, raw/REML qualification of every actual
design/control, complete timing and independently checked fits remain required.
The original covariance reader, retained qualification and timing queues are
unchanged. Ancestral outcomes/posterior assessment, accepted framework,
the atlas, inferential calibration and all eight aims remain open. GPU
inference stays paused.
