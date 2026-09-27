# Verified duplication-control selection and balance

The full metadata-based matching stage is complete and independently verified.
It produced 2,786,912 selected records across 54 sensitivity scenarios, with
873,892 target/policy records and all unmatched scenario IDs retained. Independent
enumeration checked all 47,190,168 target/policy/scenario decisions, including
eligibility, selected identities, endpoint orders, scores, ties and reuse.
Independent balance reconstruction then checked all 432 guide/policy/scenario
coverage summaries and 3,456 feature summaries.

These are completed matching diagnostics. Structural outcome qualification,
matched evolutionary effects, taxon/family weighting and phylogenetic dependence
remain unresolved; the overall project is not complete.

## Illustrative strict-background scenarios

Both examples require conservative shared architecture, native ortholog
assignments and unreported parent duplications in both guide alternatives,
sequence distance within a factor of 1.5, and moderate metadata tolerances
(endpoint length ratio ≤1.25, mean-pLDDT difference ≤10, low-confidence-fraction
difference ≤0.1 under one consistent endpoint correspondence). Counts below use
the alignment-E-value annotation policy. Equal counts across guides do not make
the guide alternatives independent replicates.

| Quantity | Any background taxon (S45) | Focal taxon required (S46) |
|---|---:|---:|
| Matched target records, each guide | 11,917 | 792 |
| Matched target taxa | 143 | 69 |
| Matched families | 2,038 | 180 |
| Distinct background nodes | 7,414 | 522 |
| Maximum reuse of a background | 22 | 22 |
| Top five backgrounds' share of matches | 0.79% | 8.59% |

The profile guide contains 109,245 original modeled targets; the mafft guide has
109,228. S45 leaves 97,328 and 97,311 unmatched, respectively. Requiring the focal
taxon further reduces coverage. These figures describe the observed modeled
terminal-duplicate universe, not all genes or all 526 sampled taxa.

## Good aggregate balance does not imply representative coverage

For profile/S45/alignment-E-value, absolute standardized target/control mean
differences range from approximately 0.00005 to 0.0885 across the eight reported
features. For example, mean pLDDT is 85.42 in selected targets and 85.52 in their
controls, and mean sequence distance is 0.298 versus 0.294. This is descriptive
balance, without an automatic pass threshold or an independent-sample test.

Selection substantially changes the target population. Mean target pLDDT is
72.25 before selection and 85.42 after selection, a shift of 0.855 original-target
SDs. Mean sequence distance falls from 0.569 to 0.298 (−0.329 SDs); length
asymmetry also decreases. Consequently, later effects in this stratum apply to
the selected, more confidently modeled conserved-architecture subset. They
cannot be generalized to all duplicates using covariate balance alone.

S45 retains 1,415 identical-model target pairs, 1,425 identical-model control
pairs and 1,447 zero-distance target/control matches. Positive-log-distance
balance uses 10,470 pairs and excludes zeros explicitly, with no epsilon.
Structural outcome processing must preserve these categories and avoid treating
identical predictions as independent evidence of biological invariance.

## Reproducibility and remaining work

Design and stage history are in [the background report](terminal-sister-backgrounds-20260927.md).
The full outputs are under:

- `results/orthology/background-control-selection-20260927-v1`
- `results/orthology/background-control-balance-20260927-v1`

Artifact hashes and full independent proofs are archived in
`metadata/background_control_selection_completed_20260927.json` and
`metadata/background_control_balance_completed_20260927.json`, with their
corresponding `_completed_readback_` files. Versioned plans and scripts reproduce
all scenarios, sparse selected rows, explicit unmatched states and diagnostics.

Next, qualify structural outcomes using verified full/domain comparisons,
coverage and confidence, retaining failures and additional selection losses.
Report contrasts with control-reuse, taxon and family dependence accounted for;
assess guide, annotation, caliper and focal-background sensitivity. Comparisons
within conserved architectures do not replace the separate domain gain/loss,
fusion/rearrangement and family-turnover analyses in the full project.
