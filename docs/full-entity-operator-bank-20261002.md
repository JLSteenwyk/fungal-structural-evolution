# Full shared-entity covariance operator bank

The complete production bank was launched on October 2 for **75,188 logical
cases**, preserving **1,052,632 endpoint/entity occurrences** from the closed
expanded covariance index. It exports **13 sparse incidence matrices** and
checks ten fixed numerical covariance benchmarks across both loading modes and
all five working species trees. Producer, independent reader and final
provenance closure are separate jobs. [Production closure is complete](../metadata/full_entity_operator_bank_completed_20261002.json);
these runs do not estimate variance components or biological effects.

The original producer finished on October 2 at 20:57 UTC. All 13 matrices,
ten full-case Grams and ten fixed numerical benchmarks are present. The
independent original reader finished at 21:25 UTC and its receipt records
`passed_full_entity_operator_raw_loading_and_gram_readback`. Final provenance
closure finished at 21:51 UTC, checking 1,830,606 source/artifact bindings and
both original journals. Both loading modes have full-bank normalized Gram
rank eight across all five working trees. The signed family kernel is zero;
the unsigned family kernel is proportional to the family contrast intercept.
Full construction plus three-right-hand-side solves took 7.95–14.70 seconds
per fixed benchmark, with maximum covariance solve residual 2.063e-12.
The original producer journal records 704.76 CPU seconds, 21.1 GiB peak memory
and no swap. Fixed benchmark timings do not measure optimization or total
project runtime.

## Explicit covariance definitions

Six namespaces retain target gene nodes, background gene nodes, physical model
pairs, families, genes and versioned model/coordinate identities. Each has the
original signed target-minus-control loadings and an unsigned reuse alternative.
Gene/model endpoints have half loadings; pair/node/family endpoints have unit
loadings. Duplicate entries coalesce and opposite signs cancel exactly.
Distinct genes remain distinct even when they share predicted coordinates.

A separate one-per-case family contrast intercept defines a family effect on
the target–control difference. It is different from signed endpoint-family
loadings, which cancel because target and control belong to the same family.
Unsigned endpoint-family incidence is twice the family-intercept incidence,
so its covariance kernel is four times the intercept kernel. These two
variance components cannot be separately identified when included together.
The bank retains both definitions and records the dependency.

### Target-node variance within the actual matching settings

The [complete original-selection proof](../metadata/full_entity_operator_target_node_redundancy_20261002_v2.json)
scanned all 4,250,692 closed selection records across 432 nonempty original
guide/policy/scenario strata. No target occurs more than once within a stratum.
Every mask, quality gate and structural screen selects a subset of those rows.
Consequently the unit target-node incidence satisfies `Z_target Z_target' = I`
within each actual setting. With a uniform residual diagonal, target-node
variance and residual variance cannot be estimated separately. The final
variance model must represent their combined contribution explicitly.

This identity does not hold across the entire pooled bank, where a target can
recur with different controls across settings. The pooled Gram rank therefore
cannot authorize a target variance within an individual cohort. Nonuniform
residual weights and other transformed loading definitions require separate
qualification. Shared genes, models and backgrounds still create dependence.
The agreement of all 51,840 producer setting-count rows was also checked;
that supplemental agreement does not replace independent input closure.
The separate full model-input stage subsequently closed with complete
Decimal/SQLite checks and original source/artifact/journal verification.

```bash
/home/bizon/anaconda3/bin/python scripts/verify_matching_target_covariance_identity.py \
  --output metadata/new-target-node-identity-proof.json
```

### Covariance identifiability after removing fixed effects

For REML, the observable variance bases are `H K_k H`, where `H = I - Q Q'`
and `Q` spans the active fixed-effect design. Raw distinct kernels can become
zero or proportional in this residual space. The new
[audit implementation](../scripts/covariance_basis_audit.py) computes both raw
and projected Frobenius Grams using component-local entity kernels and
low-rank projections. It never allocates the full covariance or `H`.
The projected Gram follows
`G_ab - 2 <K_a Q, K_b Q> + <Q' K_a Q, Q' K_b Q>`.

Zero bases, dependent bases, cancellation near numerical resolution and rank
boundaries remain explicit review states. The rounding envelope includes
contraction dimensions and normalized-design conditioning. No clipping,
jitter or automatic basis deletion changes the proposed covariance model;
deleting arbitrary dependent bases could change its nonnegative variance cone.
This is a numerical structural audit, not evidence of biological adequacy or
statistical calibration. Each production cohort and active design still needs
qualification after its source closure.

[Current software validation](../metadata/covariance_basis_audit_validation_20261002_v2.json)
passed 24 independent dense null-space comparisons, extreme design scaling,
uniform/nonuniform residual contrasts, an intercept-annihilated family kernel
and distinct raw kernels that alias under REML. It rejected 13 invalid inputs.
Maximum relative projected-Gram Frobenius discrepancy was 4.05e-16; the largest
absolute discrepancy, 3.06e-5, occurred in intentionally large-scale kernels.
An allocation check verified that the largest dense entity kernel was 6 by 6
in an 18-row, three-component test. The previous unlaunched software validation
receipt is preserved; version 2 binds the current implementation and checks.

```bash
OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python \
  scripts/check_covariance_basis_audit.py --output metadata/new-covariance-basis-validation.json
```

For supplied nonnegative variances, the working covariance is

`V = D + sum_k variance_k * Z_k Z_k' + species_variance * F F'`.

`D` is a positive residual diagonal, `Z_k` are explicit entity loadings and
`F` is the full expanded species factor for one working tree. The existing
3,330 connected family/entity components contain every shared entity; the
largest has 2,884 cases. The numerical backend factors each component and
adds the global phylogenetic term using the low-rank determinant/inverse
identities. It allocates no full 75,188-by-75,188 covariance.

The backend exposes solves, matrix application, inverse diagonals, log
determinants and direct-residual profiled REML. Solve residuals are checked;
bounded numerical refinement preserves the target covariance. Invalid
partitions, negative variances, nonfinite values, zero/under-resolved residual
quadratics and rank-deficient fixed designs are rejected. No jitter or
eigenvalue clipping repairs failures.

## Full bank and independent verification

Every original logical case remains in the bank, including cases excluded
from downstream structural masks. Cohort-specific fitting will select original
case rows through the verified design inventory. The bank neither replaces
missing responses with zero nor claims every bank case is fit-eligible.

Complete Frobenius kernel Grams compare residual identity, all six entity
kernels, the family contrast intercept and species covariance under each
mode/tree. Zero and proportional kernels remain explicit. The rank of the
normalized Gram is a diagnostic: Gram conditioning squares covariance-basis
conditioning, and qualification is still required within each actual cohort.

The producer uses sparse coordinate coalescing and component kernel inner
products. The independent reader reconstructs every raw source loading with
Counters, verifies original entity hashes and identities, and checks Gram
entries through `||Z_a' Z_b||_F^2`. Species self-products are independently
reconstructed using original pattern multiplicities. All original case order,
labels, matrix shapes and nonzero cells require exact agreement.

Ten benchmarks use fixed entity variances 0.125, a family contrast-intercept
variance 0.5, species variance 0.25 and residual diagonal one. Every case is
solved for three deterministic right-hand sides: ones, a linear ramp and a
sine of the original row index. These are numerical checks, not fitted
responses or biological estimates. The reader independently assembles each
base covariance from explicit latent outer products and uses LU for base and
phylogenetic solves/determinants. Final closure requires every source/artifact
hash and both original completion/resource journals.

## Tests and resources

The [primitive software checks](../metadata/shared_entity_covariance_fixture_validation_20261002.json)
passed 36 independent dense covariance/REML cases, row permutations, zero-rank
and zero-variance alternatives, cancelling aliases and an 80-digit strongly
correlated signed example. Maximum absolute solve, log determinant and
profiled REML discrepancies were approximately 3.4e-15, 2.2e-14 and 1.5e-14.
Thirteen invalid mathematical inputs were rejected.

The [full bank software checks](../metadata/full_entity_operator_fixture_validation_20261002.json)
passed all 13 operators and ten mode/tree benchmarks, interrupted replay and
completed-restart refusal. All 15 rehashed altered exports were rejected,
including case reorderings, collapsed genes, changed loadings, promoted zero
family covariance, wrong Grams, wrong fixed-parameter solutions and an omitted
tree. Prior source closures/journals are synthetic software fixtures, not a
biological pilot or production acceptance.

Resources were recorded before launch: **two CPU cores, 32 GiB memory, no swap
and one BLAS thread**, with 16 GiB output/scratch allowance and 100 GiB free-disk
reserve. The exact source census gives a maximum dense component of
**66,539,648 bytes**, a seven-kernel component stack of **465,777,536 bytes**,
all component-square storage of **191,079,776 bytes** and an expanded rank-301
factor of **181,052,704 bytes**. Source proof dictionaries dominate additional
memory. The original resource note's prose overestimates the first two values
by 800 and 5,600 bytes; its computed numeric fields are exact and authoritative.
The frozen note is preserved. The 1–12 hour stage planning range is
uncalibrated; measured benchmark timings will inform later fit estimates.
No GPU or paid infrastructure is used.

The [production plan](../metadata/full_entity_operator_plan_20261002.json),
[launch inventory](../metadata/full_entity_operator_bank_launches_20261002.json),
[field dictionary](../metadata/full_entity_operator_data_dictionary_20261002.tsv)
and [recorded numerical environment](../environments/shared-entity-covariance-20261002.yml)
provide source, schema and execution details. The supplementary
[launcher environment](../environments/full-entity-operator-launch-20261002.yml)
also includes the exact `psutil` version used for process identity checks.

```bash
/home/bizon/anaconda3/bin/python scripts/launch_full_entity_operators.py \
  --plan metadata/full_entity_operator_plan_20261002.json
```

This version is already launched and refuses an existing output/launch.
Reproduction requires new explicit identities while preserving original jobs.

## Remaining scientific requirements

The bank supports richer dependence calculations than reused-background and
family terms alone. It does not demonstrate that these working random effects
fully describe structural prediction uncertainty or evolutionary inheritance.
Variance identities, signed/unsigned definitions, weighting and scientific
control variants must be specified in final fit identities; unidentifiable
components must not be reported as separate estimated effects.

The [full uniform cohort/design qualification](full-uniform-covariance-qualification-20261002.md)
is now queued behind original design closure, with every setting, both loading
modes and all five trees retained. Complete cohort/design verification, expanded optimizer and numerical audits,
uncertainty/calibration and multiple-testing controls remain required. Accepted
species/reconciled gene trees, branch durations, prediction-source, domain/PAE,
missingness and ascertainment controls remain open. Working substitution
kernels are not accepted dated phylogenies. All eight scientific aims remain
incomplete; GPU prediction remains paused.
