# Domain alignment coverage sensitivity screen

This descriptive screen preserves all 70,395 domain interval pairs under both
full and pLDDT70 input masks, with both USalign input orders retained. It reads
the completed full-cohort numerical diagnostic; it does not change the original
failed strict audit or declare diagnostic output suitable for biological inference.

The analyst-defined sensitivity grid uses minimum aligned lengths of 30 or 50
residues and minimum coverage of 50%, 70% or 90%. These are exploratory screening
choices, not externally validated biological cutoffs. All six combinations are
reported. No combination was selected for a favorable duplication-effect result.

For each pair and mask, both input orders must meet the minimum aligned length
and cover the stated fraction of **each original domain interval**. Confidence
masking does not shrink the coverage denominator. Retained input lengths and
native coverage against those retained inputs are also preserved for inspection.
Threshold comparisons use exact rational arithmetic. Original coverage is
reported against interval A and interval B; native left/right metrics remain
in their original input-order orientation. No preferred order is selected and
no order metrics are averaged.

Every row retains all applicable exclusion reasons: unavailable input,
numerical discrepancy, short alignment and low original-interval coverage.
Reasons overlap, so their counts must not be summed to obtain excluded rows.
If either order has a numerical discrepancy, that pair/mask fails every screen.
Unavailable comparisons remain present with blank numeric metrics. A numerical
mismatch is retained with its original and recomputed RMSDs, never corrected
silently. Both-mask counts require the same pair to pass in full and pLDDT70
analyses, supporting later comparisons on a fixed cohort.

The independent checker compares every exported native/diagnostic metric,
interval length, retained-residue count, coverage denominator, threshold
decision, exclusion reason and aggregate against source tables. Fixtures test
exact thresholds, failure in only one order, confidence-fragment inflation,
quarantine, missing inputs, order symmetry and impossible aligned lengths.

Run from the repository root using the project Python environment:

```bash
python scripts/check_domain_coverage_threshold_cases.py
python scripts/screen_duplication_domain_alignment_coverage.py --plan metadata/duplication_domain_coverage_screen_plan_20260927.json
python scripts/check_duplication_domain_coverage_screen.py --plan metadata/duplication_domain_coverage_screen_plan_20260927.json --output metadata/duplication_domain_coverage_screen_completed_20260927.json
```

Output paths must not already exist. For a new run, use a new plan/output path;
existing source pins and artifacts are immutable.

Passing this screen establishes only numerical consistency within the diagnostic
and the stated length/coverage requirements. It does not validate domain
boundaries, establish coordinate rank, calibrate prediction error, establish
homology/function, or test structural divergence after duplication. Those
questions require the event/reference joins, sequence covariates, sampling
controls and phylogenetic analyses. Sparse or failed cases remain part of the
coverage accounting rather than becoming evidence of biological absence.

## Completed screen

All 140,790 pair/mask rows and 280,824 numeric order records passed the independent table readback. The 378 unavailable pair/masks and one numerically discrepant pair/mask remain explicit.

| Minimum aligned residues | Original coverage | Full mask | pLDDT70 mask | Both masks |
|---:|---:|---:|---:|---:|
| 30 | 50% | 70,224 | 68,270 | 68,270 |
| 30 | 70% | 69,374 | 64,036 | 64,036 |
| 30 | 90% | 61,286 | 43,390 | 43,388 |
| 50 | 50% | 65,278 | 62,747 | 62,747 |
| 50 | 70% | 64,522 | 59,059 | 59,059 |
| 50 | 90% | 56,993 | 40,144 | 40,142 |

See the [completion record](../metadata/duplication_domain_coverage_screen_completed_20260927.json) and [pinned plan](../metadata/duplication_domain_coverage_screen_plan_20260927.json).
