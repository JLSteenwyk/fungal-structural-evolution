# Robustness of common-core reference-similarity contrasts

This summary retains every oriented interval triad and all alternatives from
the independently verified common-core fits. It calculates ranges over eight
input-order combinations separately for each mask and mapping definition,
and over all 32 mask/mapping/order combinations jointly. Six length/coverage
screens give 239,580 summary rows for 7,986 oriented interval triads.

A direction label requires **every expected alternative** to pass its screen.
An incomplete group is labeled `incomplete`, even if its available values all
share a sign. The number passing, minimum, maximum and range of passing values
are retained; ranges for incomplete groups explicitly describe only that subset.
No preferred order, mask or mapping is selected and no distances are averaged.
The output preserves all triad keys for subsequent joins to event/domain links.

The signed contrast is RMSD(A, reference) − RMSD(B, reference), with identical
residue triples within each comparison. Positive means A is farther from this
reference under the particular mapping, not that A accumulated more changes
since duplication. Reference proteins are provisional and not ancestors.

All results are reported at three analyst-defined descriptive margins: 0,
0.01 and 0.1 Å. A direction is consistent beyond a margin only when the entire
range lies strictly above it or strictly below its negative. A range extending
beyond both margins changes direction; otherwise it touches or falls within
the margin band. Exact equality belongs to the band. These margins are not
prediction-error estimates, significance thresholds or validated biological
cutoffs. Tiny signs in the zero-margin summary must not be overinterpreted.

## Independently checked results

At the 30-common-residue/70%-original-interval-coverage screen, 7,025 triads
pass all 32 alternatives; 961 are incomplete.

| Descriptive margin (Å) | A farther throughout | B farther throughout | Extends beyond both directions | Touches/within margin band | Incomplete |
|---:|---:|---:|---:|---:|---:|
| 0 | 2,933 | 2,869 | 1,223 | 0 | 961 |
| 0.01 | 2,474 | 2,412 | 896 | 1,243 | 961 |
| 0.1 | 1,022 | 1,005 | 132 | 4,866 | 961 |

The full [completion record](../metadata/duplication_domain_common_core_robustness_completed_20260927.json)
verifies all 255,552 source rows and 239,580 summaries using independent
dataframe grouping. It checks complete alternative grids, all keys, counts,
extrema/ranges, completeness flags and classifications. Fixtures cover missing
alternatives, opposing signs, zero and exact-margin boundaries.

These are interval triads, not independent duplication events. Annotation
boundary/policy alternatives can share computations; guides, genes and tied
references remain linked through the common-residue event ledger. A count here
must not be substituted for event replication. The current summaries do not
establish robustness across all annotation boundaries/policies or reference
choices, and do not adjust for shared ancestry, family or sequence divergence.
Those joins and controls remain necessary for biological inference.

Reproduce using the [pinned plan](../metadata/duplication_domain_common_core_robustness_plan_20260927.json):

```bash
python scripts/check_common_core_robustness_cases.py
python scripts/summarize_common_core_robustness.py --plan metadata/duplication_domain_common_core_robustness_plan_20260927.json
python scripts/check_common_core_robustness.py --plan metadata/duplication_domain_common_core_robustness_plan_20260927.json --output metadata/duplication_domain_common_core_robustness_completed_20260927.json
```

Use new output paths for reruns; existing outputs refuse overwrite. This stage
uses existing fitted tables only and runs no GPU work or new structural fits.
