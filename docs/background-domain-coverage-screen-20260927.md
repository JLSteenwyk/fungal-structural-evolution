# Background domain coverage sensitivity

Applied the same six exploratory screens used for duplicated-domain comparisons
across all 66,929 background interval pairs. Both input orders must pass;
coverage uses the original interval length, including when low-confidence
residues have been removed. Threshold decisions use exact rational arithmetic.
No structural outcome was used to choose a favorable threshold.

| Minimum aligned residues | Original coverage | Full | pLDDT70 | Both masks |
|---:|---:|---:|---:|---:|
| 30 | 50% | 66,746 | 64,983 | 64,983 |
| 30 | 70% | 66,078 | 61,353 | 61,353 |
| 30 | 90% | 60,675 | 43,749 | 43,749 |
| 50 | 50% | 61,815 | 59,586 | 59,585 |
| 50 | 70% | 61,265 | 56,506 | 56,505 |
| 50 | 90% | 56,262 | 40,375 | 40,374 |

Full independent verification passed for all 133,858 rows and 267,246 numeric
order records. The complete table retains 133,858 pair/mask rows, including 235 unavailable
pair/masks and three pair/masks with numerical discrepancies (five directed
alignments). These exclusions can overlap other reasons. Masked comparisons
can align different residues from full inputs, so both-mask counts explicitly
require the same pair to pass both, rather than taking the smaller marginal
count. A geometry audit runs separately; passing coverage alone is insufficient
for scientific eligibility.

The producer is `scripts/screen_background_domain_coverage.py`, configured by
`metadata/background_domain_coverage_screen_plan_20260927.json`. The independent
checker `scripts/check_background_domain_coverage_screen.py` verifies every
pair/mask/order, numeric field, interval denominator, retained-position count,
threshold decision, exclusion reason and aggregate. Its final proof is
`metadata/background_domain_coverage_screen_completed_20260927.json`.
Threshold-boundary, original-denominator, quarantine, order-symmetry and
invalid-count fixtures passed for the background implementation.

These are distinct interval-pair counts, not independent evolutionary events.
The next analysis must join qualified comparisons to matched events, retain
repeated controls and boundary alternatives, and account for family and
phylogenetic dependence. No duplication effect is estimated here.
