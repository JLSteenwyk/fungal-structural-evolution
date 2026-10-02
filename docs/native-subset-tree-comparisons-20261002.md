# Actual taxon-subset tree comparisons

The complete comparison workflow is queued behind the 16 actual native
taxon-subset PMSF fits and their independent collection. It will compare every
new subset ML/consensus view with every matching projected concatenated
reference and coalescent candidate. This distinguishes sensitivity to
reinference on a retained taxon set from simply pruning an original estimate.
No production result is claimed before the original inference, collection,
comparison and independent readback finish.

## Scope and required evidence

The four cohorts were fixed before inspecting these results:

| Cohort | Taxa | Ingroup | Outgroup |
| --- | ---: | ---: | ---: |
| Exclude sparse boundary taxon | 525 | 500 | 25 |
| Exclude below 10% occupancy in both alignments | 522 | 499 | 23 |
| Exclude curated hybrids | 524 | 499 | 25 |
| Exclude hybrids and incomplete labels | 503 | 478 | 25 |

Each cohort has four actual subset fits crossing profile/MAFFT alignments and
profile/MAFFT guides, contributing eight ML/consensus views. Alongside those
are eight exactly matched projected reference views and six coalescent
alignment/support candidates: 22 views per cohort, 88 total. Every pair
involving an actual subset view is included:

| Comparison | Per cohort | Total |
| --- | ---: | ---: |
| Actual subset versus coalescent | 48 | 192 |
| Actual subset versus projected reference | 64 | 256 |
| Actual subset versus actual subset | 28 | 112 |
| All comparisons | 140 | 560 |

The [plan](../metadata/native_subset_tree_comparison_plan_20261002.json) requires
the completed [315-comparison source closure](../metadata/coalescent_tree_comparisons_completed_20261002_v2.json)
and the future complete native collection locator
`metadata/species_taxon_pmsf_collection_completed_20261001.json`.
The latter must certify all 16 runs, 32 views, 16,000 original bootstrap trees,
902,216 site profiles and both original native/collection completion journals.
The source loader verifies the entire bound archive, exact alignment/guide/
cohort identities, taxa, support tables and original tree hashes. A partial
native fit or projected reference cannot satisfy this requirement.

## Resolution, support and branch units

Every union split-presence cell, pairwise RF count, incompatible unique-split
pair and canonical four-taxon witness is retained. All 88 declared
ingroup/outgroup boundaries are measured as unrooted splits. Native consensus
trees may be unresolved; the exports retain each tree's observed split count
and the number of missing internal slots rather than fabricating resolutions.
Both distance denominators are explicit:

- `normalized_rf` is RF divided by `2 * (taxa - 3)`, the maximum for two
  resolved unrooted trees on the same cohort.
- `normalized_rf_by_observed_splits` is RF divided by the sum of observed
  internal splits. A star/star comparison has RF zero but denominator zero;
  its observed-split normalization is unavailable, not zero.

Neither normalization makes an unresolved consensus a well-supported agreement
with a resolved tree. All distances are descriptive comparisons of point
estimates on identical taxa; they are not a calibrated tree posterior.

Support screens remain separate: actual subset ML uses SH-aLRT at least 80
and empirical UFB at least 95%; actual subset consensus and projected
references use empirical UFB at least 95% with SH-aLRT unavailable; coalescent
candidates use local posterior at least 0.95. Missing metrics remain missing.
Native subset lengths are amino-acid substitutions per site, projected lengths
are original branch-path sums, and coalescent lengths are MAP coalescent
units. No numerical comparison across those units or likelihood ranking across
different site universes is performed. Overlapping conflicts are not
independent evolutionary events, and a role boundary does not establish a root.

## Independent verification and software contracts

The producer uses canonical split masks. The separate reader reparses the
actual raw trees with DendroPy, independently prunes original reference trees,
reconstructs branch metrics and support labels, and checks every RF count,
bipartition incompatibility, witness, boundary, presence cell and unavailable
denominator. Raw native ML/consensus labels are checked against the empirically
recomputed bootstrap frequencies with their reported precision. Full producer
and reader summaries must agree before closure verifies all hashes and both
original producer/reader completion journals. The handoff wrapper checks
original process/invocation identity rather than collected systemd defaults.

The full four-cohort/88-view/560-pair synthetic software contract passed.
It includes unresolved and all-star consensus trees, four genuine zero
observed-split denominators, exact/below/above support thresholds, and separate
branch units. All 15 altered exports were rejected, including replacing an
unavailable zero-denominator value with zero. The receipt is
`results/software-checks/native-subset-tree-view-grid-20261002-v1/receipt.json`,
SHA-256 `de3d92efcab78e8836a6a551e3c9ea7fde08e29dad996bb5e1de88b33ab6cecd`.
These are software checks of the full workflow; they are not a biological
pilot, completed native inference or production comparison result.

## Execution, resources and current dependency

Resources were estimated before queuing the stages: two CPU, 16 GiB RAM,
no swap, one BLAS thread, 4 GiB output allowance and 100 GiB free-disk reserve.
The conservative conflict-check bound is 560 times 522 squared. The
0.03–2 hour per-stage planning range is uncalibrated and excludes the native
inference wait. [Resource estimate](../metadata/native_subset_tree_comparison_resources_20261002.json).
No GPU or paid resource is used and existing inference jobs are unchanged.

The [runtime checkpoint](../metadata/project_runtime_checkpoint_20261002_v6.json)
verified all three original queued comparison handles and the original native
subset job. At that observation, the first 525-taxon profile guide was still
running with eight CPU threads and a 32 GiB IQ-TREE memory setting. Later
mixture stages retain their original 600 GiB setting and 650 GiB available-RAM
guard. Their wait is not included in the comparison estimate, and completion
has no guaranteed ETA.

The exact queued units are recorded in the
[producer launch](../metadata/native_subset_tree_comparison_launch_20261002.json),
[reader launch](../metadata/native_subset_tree_comparison_readback_launch_20261002.json)
and [closure launch](../metadata/native_subset_tree_comparison_closure_launch_20261002.json).
The eventual comparison output is
`results/phylogeny/full-native-subset-tree-comparisons-20261002-v1`;
its future small closure locator is
`metadata/native_subset_tree_comparisons_completed_20261002.json`.

```bash
python scripts/compare_native_subset_tree_views.py \
  --plan metadata/native_subset_tree_comparison_plan_20261002.json
python scripts/readback_native_subset_tree_views.py \
  --plan metadata/native_subset_tree_comparison_plan_20261002.json \
  --output results/phylogeny/full-native-subset-tree-comparisons-20261002-v1/readback.json
python scripts/close_native_subset_tree_views.py \
  --plan metadata/native_subset_tree_comparison_completion_plan_20261002.json
```

These commands document the already queued immutable jobs. Do not start
duplicates; reproduction uses new plans and output locations. Large tree,
bootstrap, profile and comparison archives remain outside Git.

Completion of these comparisons will establish measured sensitivity to actual
taxon-subset inference, conditional on the selected markers and models.
Gene-estimation uncertainty, marker dependence, taxon identity, model adequacy,
rooting, defensible dating and reconciliation still require qualification before
structural-rate and duplication analyses. All eight evolutionary aims remain
incomplete; GPU protein prediction remains paused.
