# Separate ancestral categorical verification

This extends the [completed scalar numerical check](independent-ancestral-scalar-readback-20261002.md)
toward the amino-acid state diagnostics needed for uncertainty-aware ancestors.
Software contracts and a complete source/resource inventory are available.
The original full production replay stopped October 2 at 22:20 EDT on an
ESS rounding boundary after 37 complete groups. A separately versioned full
replacement launched at 22:49 EDT after software and actual-source regression
checks. Its serialized reader and provenance closure are queued. No ancestral
posterior or new ancestral structure is qualified.

## Scope and method

The [complete inventory](../metadata/independent_ancestral_category_inventory_20261002.json)
retains all 405 original quartets: 403 complete, two unresolved, all 1,620 chain
identities and both original cutoffs. Production uses the 22-state alphabet
`ACDEFGHIKLMNPQRSTVWYX-`: unknown residues and gaps remain explicit. The inventory
contains 1,061,252 pattern rows at cutoff 250 and 971,274 at cutoff 500, totaling
2,032,526 pattern/cutoff rows and 44,715,572 declared indicator rows. Each cutoff
represents 22,871,876 node/anchor coordinates; repeated anchors are not independent
biological observations.

The [separate implementation](../scripts/independent_ancestral_categorical_diagnostics.py)
uses integer histograms for counts, reconstructs all declared binary state
indicators and calculates their marginal metrics with the direct-lag scalar
backend. It imports no original categorical/pattern implementation or ArviZ.
Absent states retain unresolved probability, an observed constant state retains
review, and any constant original indicator chain retains its no-metric review.
Empirical pairwise total variation between the four chains is descriptive;
it supplies neither an iid test nor posterior uncertainty.

Every indicator's counts, frequencies, disposition and defined metrics must
match the original report, at the unchanged scalar tolerance. Exactly constant
split trajectories retain explicit singular-metric review. Diagnostic agreement
is not joint alignment/state convergence. No numeric ordering of amino-acid
labels is used.

Lossless pattern grouping uses lexicographic unique full temporal rows and then
restores first-coordinate appearance order, separately from the original hash
dictionary. All chain, draw and coordinate identities remain intact. Patterns
with identical complete traces are computational reuse units, not homologous
sites, independent replicates, frequency bins or thinned trajectories.

## Software evidence

The [locked-oracle validation](../metadata/independent_ancestral_categorical_validation_20261002_v2.json)
passed 84 cases with 19/20/21/50/75/500/750 draws, including all declared labels,
unknown-only and unknown/gap variation, rare states, sticky trajectories,
alternation, chain disagreement and constant halves. All 84 simultaneous
label/code permutations preserve indicator metrics and frequencies. Reordering
the total-variation sum changed values by at most `2.23e-16`, within the explicit
`5e-16` rounding allowance. This allowance applies only to the permutation test;
same-label production descriptive counts/distances require exact equality.

Maximum defined errors against ArviZ 0.22.0 were `2.23e-16` for R-hat,
`5.01e-12` bulk ESS, `7.74e-12` tail ESS and `1.74e-17` mean MCSE. Twelve
constant-half indicator cases retain singular R-hat review. Six integer dtypes
were checked. Seven invalid state/alphabet inputs and eight altered exports
were rejected, including missing states, changed metrics, false status and an
added scientific-eligibility claim. Sixteen multidimensional pattern fixtures
reproduced every original temporal pattern, first-appearance ID and multiplicity;
three invalid pattern inputs were rejected. These are software fixtures, not a
biological pilot or complete production replay.

The [complete workflow fixture](../metadata/independent_baliphy_category_grid_validation_20261002.json)
retains all 405 synthetic groups, including 403 complete and two failed groups.
Four distinct temporal patterns use genuine locked-oracle metrics: constant,
nonconstant, unknown/gap variation and singular split trajectories. Both cutoffs
cover 3,224 pattern rows, 70,928 declared indicators and 3,224 unanchored-count
rows. All serialized rows were reconstructed, and 806 singular indicator rows
retained review. Interrupted whole-group checkpoints and compressed reports were
reused without changing bytes; completed producer and alternate-reader restarts
were refused. Twelve rehashed false exports were rejected, including omitted/
duplicate patterns, changed metrics/unknown counts, missing states, false source
digests or scientific eligibility, changed multiplicities/unanchored metrics,
hidden cutoffs, false failure dispositions and an omitted original failed group.
Its repeated fixture patterns and source journals are synthetic; this is not
the production 2,032,526-pattern result.

## Original production replay and provenance

The [35-pin full plan](../metadata/independent_baliphy_category_plan_20261002.json)
fixes every original group, cutoff, state and expected source count. The
[source loader](../scripts/independent_baliphy_category_sources.py) requires the
closed recovery overlay and complete 50,019-binding diagnostic archive, freshly
verifying the full source hashes. Each quartet array is compared with all four
selected chain arrays, their exact seeds, sample-audit/log links, node/coordinate
identities and unanchored counts. Chain arrays use the native `states` field;
quartet arrays use `values`. Auxiliary per-chain count arrays inherit their prior
closed projection checks, rather than a newly independent native parser claim.

The [replay engine](../scripts/independent_baliphy_category_replay.py) independently
rebuilds every retained temporal pattern, first-appearance coordinate map and
multiplicity. Every original pattern report is compared state by state,
including all absent states and the six descriptive chain distances. All four
unanchored-residue count screens per cutoff are separately recalculated, totaling
3,224 additional rows. Comparisons retain the original scalar thresholds and
`1e-8` absolute-plus-relative tolerance.

The [producer](../scripts/prepare_independent_baliphy_category_readback.py) writes
deterministic gzip JSONL with the complete separately calculated indicator report,
original diagnostic digest, defined errors and unresolved flags for every
pattern. Cutoff files publish atomically after flush/fsync; whole-quartet JSON
checkpoints bind both cutoffs. Interrupted files resume by recalculating every
saved row without rewriting accepted files. Completed producers refuse restart.

The [serialized reader](../scripts/readback_independent_baliphy_category_readback.py)
rebuilds every array/pattern/indicator comparison and exact serialized row,
reconciles all source partitions, rejects foreign/missing/extra artifacts and
rechecks full source/output hashes. It uses the same new estimator and is not a
third independent numerical implementation. Final closure requires both actual
original producer/readback completion journals and the full source/artifact graph.

The [three-stage launch inventory](../metadata/independent_baliphy_category_launches_20261002.json)
records the original producer controller 2670368, reader 2670372 and closure
2670376. The first [execution observation](../metadata/independent_baliphy_category_execution_checkpoint_20261002.json)
confirmed all three live original identities/cgroup limits and three completed
producer checkpoints: 12,452 pattern rows and 273,944 indicators. These are
unclosed producer progress, not final verification or posterior qualification.

## Sources, resources and remaining work

The [inventory script](../scripts/inventory_independent_ancestral_categories.py)
rechecks small manifests, report receipts, summaries and the closed diagnostic
archive hash. All group identities, four distinct seeds, 101 original saved
iterations, 75/50 retained draws, 22 states and coordinate/status totals reconcile
with the existing 50,019-binding closure. Large data hashes are copied from that
closed archive and are explicitly **not rehashed or numerically replayed** by
the inventory. Complete launch/readback verification must check them afresh.

Observed compressed pattern arrays occupy 122,980,341 bytes, categorical JSONL
reports 446,391,252 bytes, and uncompressed full quartet state arrays total
9,240,237,904 bytes. The largest quartet array is 460,547,072 bytes and has
1,139,968 coordinates; the largest retained pattern bank has 59,810 patterns.
An estimated sixteen-array workspace is 7,368,753,152 bytes, not a hard bound.
Installed per-stage limits are two CPU equivalents, 32 GiB, no swap, one BLAS thread,
128 GiB output allowance and 228 GiB minimum free disk. The 4–96 hour range per
stage is uncalibrated planning, not measured runtime or a project finish ETA.
No sampling, GPU inference or paid resource was launched by this inventory.

The original full replay, serialized readback and actual-invocation/hash closure
must finish before a complete production numerical-validation claim. Independent native
alignment parsing, adequate sampling horizons, joint posterior/model/root
qualification and predictor controls remain separate requirements. Both native
allocation failures remain unresolved; GPU inference is paused. Aim 8 and all
other biological aims remain incomplete.

Reproduce software contracts in the locked ancestral diagnostics environment:

```bash
SOFTWARE/ancestral-diagnostics-20260927/bin/python \
  scripts/check_independent_ancestral_categorical_diagnostics.py --output NEW.json
/home/bizon/anaconda3/bin/python scripts/inventory_independent_ancestral_categories.py \
  --output NEW_INVENTORY.json
SOFTWARE/ancestral-diagnostics-20260927/bin/python \
  scripts/check_independent_baliphy_category_readback.py --output NEW_GRID_CHECK.json
```

Production scripts require `--plan` and the reader also takes `--output`.
Reproduction requires new explicit plan/output/invocation identities; do not
edit or relaunch the frozen production plan. The reproduction environment remains
[Python 3.10.13/NumPy 2.2.6/SciPy 1.15.3](../environments/independent-ancestral-verification-20261002.yml),
with the separate locked ArviZ environment used only for oracle fixtures.


## Exact ESS boundary and full replacement, October 2 evening

The [three original failure journals](../metadata/independent_baliphy_category_v1_failure_20261002.json)
remain preserved with the original 35-pin plan, code, 37 complete group
checkpoints and partial next group. Producer controller 2670368 failed at
22:20:19 EDT; reader 2670372 and closure 2670376 failed on that dependency.
This differs from the earlier scalar tied-percentile problem and does not
indicate a new native allocation failure.

The [actual-source regression check](../metadata/independent_baliphy_category_actual_regression_check_20261002.json)
freshly hashed the original categorical pattern bank and diagnostic gzip,
extracted pattern 932 at cutoff 500, and reproduced the old assertion:
K-state bulk ESS 188.67018722995684 versus 199.59370238699853. The R-hat,
tail ESS, mean MCSE and original mixing disposition agree. A small exact trace
is [versioned for reproducible software checks](../metadata/binary_ess_boundary_regression_20261002.json).

For this binary indicator, exact integer/Fraction arithmetic gives lag-2 and
lag-3 autocorrelations 38/655 and −38/655. Their paired sum is exactly zero.
The FFT path rounds the pair to zero while direct arithmetic rounds it to a
tiny positive value. Strict positive-pair truncation then stops at different
lags. The [certificate implementation](../scripts/binary_ess_boundary_certificate.py)
enumerates stopping/storage alternatives only at reachable exactly zero
pairs, terminating at a genuinely negative pair. Both observed ESS values
must match an admissible exact outcome at the unchanged 1e-8 absolute-plus-
relative tolerance. Other mismatches still fail. Differing certified metrics
remain `exact_binary_zero_pair_truncation_requires_review`, are excluded from
ordinary-agreement error maxima, and are counted separately from singular
R-hat reviews. This is unresolved numerical review, not numerical agreement
or posterior acceptance.

The [certificate tests](../metadata/binary_ess_boundary_validation_20261002_v2.json)
passed the real failing trace and 250 random binary fixtures; +1 alterations
to R-hat, bulk ESS, tail ESS and MCSE were rejected. The [replacement full-grid
workflow test](../metadata/independent_baliphy_category_v2_grid_validation_20261002.json)
passed all 405 synthetic groups, retaining two failures, 4,030 pattern rows,
88,660 declared indicators and 3,224 unanchored screens. Five shared patterns
include the exact real regression: 403 ESS-boundary indicator reviews and 806
singular reviews remain explicit. All serialized reconstruction, byte-identical
interrupted-checkpoint reuse, completed-restart refusal and 12 rehashed
false-export rejection checks passed. These repeated software fixtures are
not the complete production result.

The [54-pin replacement plan](../metadata/independent_baliphy_category_v2_plan_20261002.json)
retains all 405 original groups, both cutoffs, all 22 states, all 2,032,526
pattern rows and all 44,715,572 declared indicators. Separate output root
`results/ancestral/full-independent-categorical-readback-20261002-v2` preserves
the failed v1 output. Its [launch inventory](../metadata/independent_baliphy_category_v2_launches_20261002.json)
records producer controller 2701791, reader 2701795 and closure 2701799.
[First execution observation](../metadata/independent_baliphy_category_v2_execution_checkpoint_20261002.json)
verified the original live handles, native producer and installed two-CPU/
32-GiB/no-swap limits. No group checkpoint was yet complete at that observation;
full production replay, serialized reconstruction and provenance closure remain
pending. At 22:53 EDT the [next observation](../metadata/independent_baliphy_category_v2_execution_checkpoint_20261002_v2.json)
confirmed four complete producer checkpoints, 28,021 pattern rows and 616,462
indicators; these are unclosed progress and do not certify complete production.
The 128-GiB output allowance and uncalibrated 4–96-hour range per stage
are planning assumptions, not a finish ETA.

Reproduce the replacement checks and inspect its frozen plan with:

```bash
SOFTWARE/ancestral-diagnostics-20260927/bin/python \
  scripts/check_binary_ess_boundary_certificate.py --output NEW_BOUNDARY_CHECK.json
SOFTWARE/ancestral-diagnostics-20260927/bin/python \
  scripts/check_independent_baliphy_category_readback_v2.py --output NEW_GRID_CHECK.json
python scripts/record_independent_baliphy_category_checkpoint.py \
  --plan metadata/independent_baliphy_category_v2_plan_20261002.json \
  --output NEW_EXECUTION_OBSERVATION.json
```

The existing three production invocations must not be relaunched. Verification
is CPU-only; native sampling and GPU prediction were not restarted.

## Independent native decoding inventory

A [separate manual FASTA/Newick/alignment projection implementation](../scripts/independent_native_ancestral_alignment.py)
imports neither Biopython nor the original native-block parser/anchor builder.
Its [software checks](../metadata/independent_native_ancestral_alignment_validation_20261002.json)
passed 80 randomized projections, explicit unknown/unanchored fixtures, five
Newick grammar cases, a 2,399-node iterative tree and native iteration fixtures.
Twenty-three malformed inputs were rejected. Expected input X residues retain
wildcard semantics; gaps, duplicated tip anchors and unanchored candidate
residues remain explicit.

The [full source inventory](../metadata/independent_native_alignment_inventory_20261002.json)
completed under original controller 2685636/native process 2685644 with
[actual completion-journal evidence](../metadata/independent_native_alignment_inventory_completed_20261002.json).
It retained all 1,620 identities, checked 1,618 intact chains including the six
intact chains belonging to unresolved quartets, and preserved the two failures.
Every original input FASTA and all runtime tree labels/tips were independently
decoded; all 50,066 source bindings were freshly hashed. Detailed source records
and the hash archive remain outside Git with paths and digests in the compact
inventory.

This census covers 45,530,706,578 bytes of native alignment files and identifies
163,418 saved alignments, 653,672 candidate frames, 9,732,673,504 state observations
and both cutoff count arrays for future independent raw-sample replay. It did
not independently decode those saved blocks or recompute their projected states,
unanchored counts or count arrays. Source-to-runtime clade mappings inherit the
closed audit; this is not independent topology or posterior qualification.
The proposed complete replay uses two CPUs/32 GiB/no swap, a 4-GiB output
allowance and 104-GiB minimum free disk; its 1–48-hour range per stage is
uncalibrated and no replay has launched. Adequate sampling, accepted phylogenetic/
model/root/predictor controls and all eight biological aims remain incomplete.


## Actual formerly failing group passed, October 2 at 23:12 EDT

The [fresh actual-source/output checkpoint](../metadata/independent_baliphy_category_v2_boundary_checkpoint_20261002.json)
rehashed the original pattern/report and the new closed group/cutoff artifact.
It independently reconstructed pattern 932 in the originally failing group
and matched the exact saved row. Both differing ESS values and their exact
zero-pair certificate remain explicit numerical review; this is not numerical
agreement. At that observation the new producer had 48 complete group
checkpoints/325,551 pattern rows/7,162,122 declared indicators. Full serialized
readback/provenance closure remains pending.

The [complete independent native sample replay](independent-native-alignment-replay-20261002.md)
has now launched after its full 1,620-chain workflow contracts passed. This
supersedes the inventory's historical unlaunched proposal. It preserves all
1,618 intact chains, including six in unresolved quartets, and both failures,
and independently checks every saved alignment/state/count cell. Adequate
joint sampling and biological posterior/model/root/predictor qualification
remain open.
