# Separate ancestral categorical verification

This extends the [completed scalar numerical check](independent-ancestral-scalar-readback-20261002.md)
toward the amino-acid state diagnostics needed for uncertainty-aware ancestors.
Software contracts and a complete source/resource inventory are available.
Production categorical numerical replay has not started, and no ancestral
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
Proposed resources are two CPU equivalents, 32 GiB, no swap, one BLAS thread,
128 GiB output allowance and 228 GiB minimum free disk. The 4–96 hour range per
stage is uncalibrated planning, not measured runtime or a project finish ETA.
No sampling, GPU inference or paid resource was launched by this inventory.

Next, implement immutable complete-grid checkpoints, full raw-array/coordinate
reconstruction and every pattern/state numerical comparison, followed by full
serialized readback and actual-invocation/hash closure. Independent native
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
```
