# Complete expanded design and recipe inventory

The full inventory was queued on October 2 behind the original model-input
closure. It covers **622,080 fixed model-setting records** and five working
species-tree alternatives, giving **3,110,400 nominal tree-setting inputs**.
These are input counts before exact sharing and before final dependence,
weighting and control variants; they are not completed fits or a final count
of all scientific models.

The [production plan](../metadata/full_expanded_model_designs_plan_20261002_v2.json)
and [launch inventory](../metadata/full_expanded_model_designs_v2_launches_20261002.json)
record immutable source specifications and the original producer, independent
reader and closure handles. Full production acceptance remains pending. A rank-boundary test exposed the
original checker's fixed relative condition-number comparison. The three
original queued invocations were [stopped and journal-preserved before production](../metadata/full_expanded_model_designs_superseded_attempt_20261002_v1.json);
new version 2 passed the corrected boundary contract before launch. Original
source files, plans and launch records remain unchanged.

## Complete grid and exact sharing

Every setting from the complete input stage is expanded across two responses,
two sequence axes and polynomial degrees one, two and three. Both guides, four
matching policies, all 54 scenarios, both residue masks, both eligibility
gates, six original screens and five physical order contrasts remain mapped.
Empty and unestimable settings remain in the table.

A cohort consists of original logical cases sorted by their full case IDs,
within one guide and residue mask. Multiple original settings can share a
cohort only when these ordered case IDs agree exactly. Distinct genes sharing
physical predictions remain distinct cases. A design further specifies order,
sequence axis, degree, predictor definitions and the closed input/covariance
source contract. A response-input identity also specifies the ordered observed
response values. Different masks, orders, scientific definitions or source
contracts are never merged merely because some numerical values happen to
match.

Each cohort stores original case row indices once. Designs and response
inputs reference the original ten Parquet partitions, avoiding repeated dense
exports. The complete setting table maps every requested specification to
these exact inputs. Shared computation does not create independent
replication or justify treating sensitivity settings as separate samples.

## Declared terms, scaling and rank

Each design declares an intercept, the selected first through third sequence
power contrasts, and all six coverage/length/confidence contrasts. Powers were
computed before averaging physical orders in the closed input stage. Columns
are not centered, preserving the zero-contrast interpretation of the
intercept.

For stable numerical rank diagnostics, each column is divided first by its
maximum absolute value and then by its L2 norm. Scaling is recorded in two
parts to avoid underflow or overflow from multiplying them. A declared column
whose values are exactly zero is retained in the audit but inactive in the
estimable matrix. No approximate cutoff or deletion of nonzero terms is used
to repair rank. In particular, the aligned pLDDT-70 fraction contrast should
be zero under the pLDDT-70 mask, where all aligned residue pairs meet the
threshold on both sides.

The producer computes singular values of the active unit-L2 design. Its
default numerical tolerance is `max(N, nominal columns) * float64 epsilon *
largest singular value`. Ranks at one tenth and ten times that tolerance are
also recorded. A disagreement between these ranks requires review. An
independent reader recomputes normalization in long-double arithmetic and uses
pivoted QR followed by LAPACK `gesvd` to check the singular values and ranks.
The active indices and all exact-zero indices remain explicit. The exported
condition number is checked against its saved singular-value ratio; its
reciprocal is independently checked within a dimension/epsilon numerical
error bound. This avoids requiring an unjustified fixed relative agreement
between ill-conditioned condition numbers while leaving singular-value checks,
rank thresholds and review dispositions unchanged.

| Disposition | Interpretation |
|---|---|
| `empty_setting` | No eligible observations; retained in the full grid |
| `sequence_axis_uninformative_requires_review` | Every declared sequence contrast is exactly zero |
| `insufficient_residual_dimension` | Observations do not exceed active coefficients |
| `rank_boundary_requires_review` | Numerical ranks depend on the recorded tolerance band |
| `rank_deficient_requires_review` | Nonzero columns remain linearly dependent |
| `constant_response_requires_review` | The observed response is constant despite a full-rank design |
| `ready_for_working_covariance_fit` | Input rank/response checks pass; model adequacy and inference are not established |

The stage does not add jitter, change thresholds after inspecting outcomes,
impute missing measurements or promote review states to accepted results.

## Readback, tests and resources

SQLite independently reconstructs each cohort from all 4,250,692 original
selection links, enforcing one selected control per target within each
original setting. Cohort order and membership are checked against every
saved NPZ. The reader verifies all raw/scaled matrix and response hashes,
source/covariance contracts, background gene-node and physical-pair reuse,
family component and species-pattern identities, and all complete setting
foreign keys. Final closure requires complete source/artifact checks and both
original process completion/resource journals.

The [software receipt](../metadata/full_expanded_model_design_fixture_validation_20261002_v2.json)
records full positive contracts for 48 distinct gene-context cases and 3,840
fixed model settings, with empty, uninformative, underspecified, rank-deficient and constant-response
states in the full synthetic grid and a separate numerical rank-boundary
contract exercised. It verifies shared-cohort computation,
full grid retention, scales as extreme as `1e-200` and `1e200`, interrupted
replay and completed-restart refusal. All **23 rehashed altered exports** were
rejected, including changed ranks/scales, removed settings, collapsed or
reordered cases, broken dependence hashes, changed trees and promotion of
constant responses. Prior closure/journal fixtures are synthetic; these tests
are neither a biological pilot nor production acceptance.

Before launch, resources were fixed at **two CPU cores, 32 GiB memory, no swap
and one BLAS thread**, with a 32 GiB output/scratch allowance and 100 GiB disk
reserve. The completed input producer bounds all cohort-member occurrences
before sharing at **60,014,752**, or at most **480,118,016 bytes** of int64
indices. The producer's 300,073,760 retained setting occurrences include five
correlated order contrasts and must not be read as effective sample size.
The inventory streams matrices with at most ten columns and allocates no
observation-by-observation covariance. Full source maps and SQLite verification
dominate additional memory/disk. The 1–24 hour planning range per stage is
uncalibrated. No GPU or paid infrastructure is used.

```bash
/home/bizon/anaconda3/bin/python scripts/launch_full_expanded_model_designs_v2.py \
  --plan metadata/full_expanded_model_designs_plan_20261002_v2.json
```

This frozen version is already queued. Reproduction requires new explicit
output and launch identities, preserving original jobs and evidence.
The [field dictionary](../metadata/full_expanded_model_designs_data_dictionary_20261002.tsv)
defines all exported fields.

## Remaining scientific work

The original producer and independent reader have both completed the full
622,080-setting inventory: 4,340 cohorts, 130,200 designs and 260,400 response
inputs, with every reported design full rank and every fit input source-ready.
The reader's original completion journal records 20:30 EDT on October 2.
The original closure handle is still live and checking the complete source/
artifact graph; these counts do not yet establish a closed production handoff.
The [queue checkpoint](../metadata/full_shared_entity_timing_queue_checkpoint_20261002.json)
verifies both exact terminal journals and the original live closure, with
qualification and full-scope timing still waiting on their original dependencies.
No production fit or calibrated inference follows from these rank results.

Once full readback closes, actual unique input counts and rank dispositions
will inform the complete covariance-fit resource estimate. A fit-input ID
identifies data and fixed-effect design; final model identities must also bind
tree, signed or unsigned shared-entity loadings, variance structure, weighting
and other prespecified controls. The shared family/entity incidence and five
working rank-301 species factors are already source-bound but are not fitted
here.

Numerical estimability does not establish ancestry correction, residual
adequacy, calibrated uncertainty or a biological association. Accepted species
and reconciled gene phylogenies, branch-duration uncertainty, prediction-source,
domain/PAE, missingness/uneven sampling/ascertainment and multiple-testing
controls remain required. All eight scientific aims remain incomplete.
