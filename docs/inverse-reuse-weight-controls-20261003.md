# Original-cohort inverse-reuse sensitivity controls

The full exporter and independent reader are software qualified and launched
for all **4,340 original cohorts**, 75,188 logical cases and 34,110,120
cohort-row occurrences. Four controls give 136,440,480 case/control
occurrences. The real exporter has written every cohort and exited
successfully; independent SQL/Fraction readback is running and complete
two-journal closure remains pending. These are prepared working-model inputs,
not fitted evolutionary results.

## Definition and interpretation

The controls are uniform, inverse background-node reuse, inverse physical
background-pair reuse and inverse connected-family-component reuse. The
group keys come from the closed original case/covariance tables. Physical
pairs use the original versioned `background_pair_key`; connected components
use `family_component`, rather than treating each family identifier as an
independent block. Counts are reconstructed within each original cohort after
its existing matching and quality gates. A logical case contributes once;
`selection_records` and setting multiplicity do not expand the observations.

For a cohort of size `n`, let `G` be its number of represented groups and
`c_i` the number of cohort members sharing observation `i`'s group. Define:

```text
w_i = n / (G * c_i)
D_i = (G * c_i) / n
```

The exact weight sum is `n`, the mean weight is one, and each group receives
total weight `n/G`. The reciprocal diagonal generally does not have mean one.
Uniform baseline weights and diagonals are one. Integer counts, group totals,
normalization numerators/denominators and case-row identities remain available
so every floating-point value can be independently reconstructed.

These are prespecified reuse sensitivity controls. Using `D` as a residual
variance kernel is a working covariance assumption; it is not a statement
that literal weighted likelihoods or balanced losses have the same semantics.
The controls are neither inverse-probability weights nor calibrated prediction
precision. Group counts and weight sums are not effective sample sizes,
independent biological observations or posterior ESS. No confidence-to-error-
variance transformation is used.

## Covariance and source requirements

When `G*c_i == n` for every row, the control is exactly uniform-one; that
classification is based on integer identities, not a numerical rank test.
Uniform controls require the named uniform fold and fresh numerical
qualification. Nonconstant controls require the already closed
[positive-diagonal cone](nonuniform-covariance-cones-20261003.md), retaining
residual `D` and target identity `I` separately, plus fresh qualification.
This exporter establishes neither weighted raw/REML identifiability nor
conditioning, whitening or precision envelopes. Near-uniform controls may
remain numerically unresolved. No term is deleted on a rank result.

Source checks verify the five closed parent completion/archive identities,
then freshly hash every consumed original case table, covariance table,
case-ID order, cohort manifest/membership NPZ and cone membership record.
They check case-table joins and original order, both cone modes for every
cohort, guide/mask identities, sorted local IDs and row/hash/count consistency.
This is explicit scoped verification: unrelated archived design/response
arrays, operators and native outputs are not freshly replayed here. Their
original broader closures remain prerequisites for numerical fitting.

Each cohort exports its original row vector, three integer reuse-count
vectors, four weight vectors and four reciprocal-diagonal vectors. Cohort
receipts bind every array digest, source contract, policy definition and
original membership. The manifest binds all cohort receipts and NPZ files.
The independent reader uses SQL partition counts and exact Fraction division,
checks every saved numerical cell and metadata field, reconstructs the full
policy census and verifies all consumed sources and output hashes. Completion
requires both actual original process journals and complete stage provenance.

## Qualification and resources

Software qualification passed 40 independent grouping cases including single
groups, unique groups, equal reuse and unequal/random reuse. It covers the
complete declared 24-case, five-cohort synthetic source/export/readback,
including uniform and nonuniform classifications for all three reuse policies.
It rejected 17 rehashed output alterations, six rehashed source alterations
and three invalid primitive inputs. Completed producer/reader restarts were
refused; positive outputs and source bytes were restored exactly and rehashed.
Parent closure and journal fixtures are explicitly synthetic. These are
software checks, not a biological pilot or full real-source result.

The real producer census reports all three reuse controls as nonuniform in
every cohort. Background-node and physical-pair maximum reuse is 71;
connected-family-component maximum reuse is 771. The background controls
have weights ranging from approximately 0.02085 to 1.5964, and the component
control from 0.00795 to 7.4521 across all cohorts. These are preliminary
export descriptions pending independent raw-count readback, not fitted
effects or an estimate of independent biological observations.

A complete post-export audit freshly verified all 4,384 consumed source
bindings and 8,681 artifacts with original producer terminal-success evidence.
It compared background-node and physical-pair weights/diagonals over every
exported row: the two controls are exactly equal in all 4,340 cohorts. Their
distinct policies and provenance remain retained. Once independent raw-count
readback passes, future numerical computation can consider exact input reuse
without treating duplicate controls as independent evidence. This byte and
redundancy audit does not independently reconstruct raw source counts or
close the full control stage. Its actual original audit wait and exact
PID/create/command/invocation journals exited zero.

The actual original software wait exited zero, with matching exact wrapper
PID/create/command/invocation journal messages and native receipt hashes.
The checker used 1.09 child CPU seconds and reported peak child RSS of
103,809,024 bytes. Manager memory accounting is distinct from native RSS.

The full stage has two CPUs, 16 GiB memory, no swap and one BLAS thread.
Each native stage has a 12-GiB address-space limit, 21,600-CPU-second budget
and 256-MiB per-file limit. Raw arrays budget 3,274,571,520 bytes (about
3.05 GiB), with at most 2,201,664 raw bytes per cohort. Metadata budgets
142,213,120 bytes; output allowance is 8 GiB and minimum free disk 128 GiB.
The provisional 60–14,400-second range per stage is uncalibrated. Budgets
and output allowances are not runtime measurements or directory quotas.
No GPU prediction, paid service or failed original job is restarted.

## Evidence and reproduction

- [Software validation](../metadata/inverse_reuse_weight_software_validation_20261003_v1.json)
  and [actual original wait/journal proof](../metadata/inverse_reuse_weight_software_transport_20261003_v1.json).
- [Full immutable plan](../metadata/full_inverse_reuse_weight_plan_20261003_v1.json),
  [prelaunch resources](../metadata/full_inverse_reuse_weight_resources_20261003_v1.json)
  and [original launch inventory](../metadata/full_inverse_reuse_weight_launches_20261003_v1.json).
- [Current original runtime](../metadata/full_inverse_reuse_weight_execution_checkpoint_20261003_v1_2256.json)
  and [field/array dictionary](../metadata/inverse_reuse_weight_data_dictionary_20261003.tsv).
- [Complete export-byte and control-redundancy audit](../metadata/full_inverse_reuse_weight_export_verification_20261003_v1.json)
  and [actual original audit wait/journal proof](../metadata/full_inverse_reuse_weight_export_audit_transport_20261003_v1.json).

`check_full_inverse_reuse_weights.py` reproduces the complete declared
software qualification under the saved software resource/wrapper configuration
with fresh output and receipt paths. `run_full_inverse_reuse_weights.py`
has separate producer and `--reader` modes. Restore the immutable pinned
versions and source archives; use fresh plan/output paths on another machine.
Completed roots are refused. Incomplete same-plan export checkpoints are
checked before reuse, but failed original scientific attempts are preserved
and are not restarted by this stage. The saved original commands record the
enforced environment and resources.

Full independent real readback/closure, weighted raw/REML numerical qualification,
complete timing and independently checked fits remain required. Inferential
calibration, accepted phylogeny/reconciliation, adequate ancestral posteriors,
the complete atlas and all eight scientific aims remain open. GPU inference
stays paused.
