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

## Balance and selection figure

![Matched balance and selection shifts](figures/background_control_balance_20260927.png)

Panel A shows selected-target minus control means divided by pooled matched SD;
panel B shows selected-target minus original-target means divided by the original
target SD. The axes share a range but use different denominators. Both guides
are sensitivity alternatives, not independent replicates. S45 is illustrative,
not a designated confirmatory primary analysis. Positive-log distance excludes
zero-distance matches. These are descriptive diagnostics, not structural effects.

Download [PDF](figures/background_control_balance_20260927.pdf),
[SVG](figures/background_control_balance_20260927.svg), or
[plotted values](figures/background_control_balance_20260927.tsv). All 32 values
and sample counts were checked against the fully audited source table; pooled
matched SD differences were additionally reconstructed from moments. Export
hashes, reproduction commands and visual review are recorded in
`metadata/background_control_balance_figure_*_20260927.json`.

## Taxon and family concentration

Full independent reconstruction passed all 2,786,912 selected records, 513,240
represented-group rows and 864 summaries (432 configurations × two grouping
levels). Group tables include represented groups, their original target counts
and selection fractions. Summaries explicitly count original groups with no
matches; the original population is all modeled terminal targets per guide.

For profile/alignment-E-value:

| Diagnostic | S45: any background taxon | S46: focal taxon required |
|---|---:|---:|
| Top five taxa, share of selected targets | 35.62% | 36.87% |
| Top five families, share of selected targets | 12.16% | 33.84% |
| Largest taxon contribution | 1,595 | 83 |
| Largest family contribution | 593 | 142 |
| Original taxa with no selected targets (of 153) | 10 | 84 |
| Original families with no selected targets (of 20,492) | 18,454 | 20,312 |

S45's 143 represented taxa have an inverse concentration of 25.05, calculated as
N² / sum(group count²). This describes the equivalent number of equally
represented groups; it is **not** an independent sample size or a phylogenetic
correction. Its 2,038 families have an analogous value of 170.33. These results
motivate taxon/family weighting and sensitivity checks alongside explicit
phylogenetic modeling; confidence balance alone does not address dependence.
Weighting changes the population-average contrast and must be stated when
reporting effects. No structural outcomes were used in these diagnostics.

Full outputs: `results/orthology/selected-control-concentration-20260927-v1`.
Plans, script hashes, reproduction commands and independent proof are archived
in `metadata/selected_control_concentration_*_20260927.json`.

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

## Structural coverage for the original matching targets

All 218,473 modeled target nodes (109,245 profile-guide and 109,228 MAFFT-guide)
now have the same six original-protein coverage screens applied to their audited
structural measurements as the background candidates. All 436,946 target/mask
rows and 2,621,676 decisions were independently replayed with decimal ceiling
cutoffs. Both native orders must pass; confidence-mask coverage uses the original
protein length. The 6,354 identical-model targets per guide remain explicitly
unmeasured rather than being assigned zero structural divergence.

Requiring 50 residues and both masks retains 43,515/43,517 profile/MAFFT targets
at 50% coverage, 26,192/26,196 at 70%, and 7,771/7,772 at 90%. These counts refer
to all original modeled targets, not only the selected matches; guides are
sensitivity alternatives rather than independent replicates. The original
control choices remain unchanged. The next step is to measure joint target/control
attrition within every fixed guide/policy/scenario before estimating effects.

Reproduce with `scripts/screen_matched_duplication_target_coverage.py` then
`scripts/readback_matched_duplication_target_coverage.py`. Output:
`results/structural_comparisons/matched-target-coverage-20260928-v1`.
Evidence: `metadata/matched_target_coverage_completed_20260928.json`.

## Full fixed-match structural attrition running

`fungal-fixed-match-structural-attrition-20260928.service` applies the completed
target and background coverage results to all 2,786,912 frozen selections.
Original selections are preserved verbatim with added target/background pass
bitsets; bit positions are bound to the six recorded screen definitions.
Full, confidence-masked, and both-mask intersection results are retained.
All 432 guide/policy/scenario strata will have 18 cells (7,776 total), recording
original targets, metadata matched/unmatched counts, and both-pass,
target-only-pass, background-only-pass and neither-pass counts. Empty strata
remain explicit. Background matching-node identities are joined to the exact
original gene pair and verified against structural pair identities.

No controls are reselected after seeing structural quality. Attrition can change
the target population; these remain dependent, descriptively matched records.
Full output and independent readback are still pending. Resource estimate:
one CPU, 16 GiB RAM, no swap/GPU or charges, up to 2 GiB output, 0.1–2 hours.
Script: `scripts/assess_fixed_match_structural_attrition.py`; execution plan and
exact process identity are in `metadata/fixed_match_structural_attrition_plan_20260928.json`
and `metadata/fixed_match_structural_attrition_launch_20260928.json`.

The fixed-match attrition producer completed successfully. A full independent
readback is now running as `fungal-fixed-match-structural-attrition-readback-20260928.service`.
It compares every original selection field, reconstructs all six target/control
bitsets from the original coverage sources and exact matching-node identities,
and independently accumulates 18 × 4 outcome-category arrays per matching
stratum. All 7,776 summary cells, including empty strata and metadata-unmatched
denominators, must agree. The checker imports no producer selection/aggregation
helpers. Production completion is not yet accepted as a validated attrition result.
Resources: one CPU, 16 GiB RAM, no swap/GPU/charges, 0.1–2 hours and negligible
proof output. Script and launch provenance are recorded under
`readback_fixed_match_structural_attrition.py` and the matching September 28
readback plan/launch metadata.

## Fixed-match structural attrition independently verified

The producer and full readback both exited successfully. All 2,786,912 original
selection records, 16,721,472 eligibility bitsets, and 7,776 summary cells passed.
Completion evidence is in `metadata/fixed_match_structural_attrition_completed_20260928.json`
and `metadata/fixed_match_structural_attrition_completed_readback_20260928.json`.

For the illustrative alignment-E-value strict-background scenarios above,
requiring both masks and both orders to align at least 50 residues covering 70%
of both original proteins gives identical descriptive counts in the two guide
alternatives (which are not independent replicates):

| Structural eligibility | S45: any background taxon | S46: focal taxon required |
|---|---:|---:|
| Original metadata-matched records | 11,917 | 792 |
| Target and control both pass | 7,914 | 656 |
| Only target passes | 387 | 21 |
| Only control passes | 588 | 37 |
| Neither passes | 3,028 | 78 |

Thus 66.4% of S45 matches and 82.8% of S46 matches survive this particular
screen. The original modeled-target denominators remain 109,245/109,228 by
guide; metadata-unmatched records remain explicit. No replacement control was
selected after screening. Surviving matches need fresh balance/representation
assessment and phylogenetic/family-aware effect estimation; these counts alone
do not establish a duplication effect or represent all sampled taxa.

## Balance after structural screening running

A full assessment now covers all 7,776 retained scenario/mask/screen cells,
including empty cells, with 62,208 rows across the same eight covariates.
Besides target/control balance and shifts from the original modeled-target
pool, it records the additional shift relative to the pre-screen metadata-matched
target pool. Retained taxa/families, distinct controls, maximum reuse, top-five
reuse share and zero sequence distances are retained. No control is reselected.

The established balance helper passed known moments/SMD/selection-shift fixtures,
empty/single/zero-variance cases, positive-log exclusions and endpoint swaps.
Full production and independent output readback remain pending. Resources: one
CPU, 24 GiB RAM, no swap/GPU/charges, up to 1 GiB output and a broad 0.1–8-hour
estimate. The executable stage is `scripts/assess_screened_match_balance.py`;
plan and exact live identity are in `metadata/screened_match_balance_plan_20260928.json`
and `metadata/screened_match_balance_launch_20260928.json`.

### Producer complete; independent reconstruction running

The producer finished successfully September 28 at 22:34 EDT, with all
2,786,912 original selections represented in 7,776 coverage cells and 62,208
feature summaries. Source and output hashes were verified after terminal exit
zero; the receipt is recorded in
`metadata/screened_match_balance_producer_completed_20260928.json`.

`scripts/readback_screened_match_balance.py` now independently reconstructs every
cell from the audited fixed-match membership and original graph nodes. It uses
separate feature/statistics implementations, decodes each endpoint's eligibility
before intersection, and checks both target baselines, moments, quantiles,
nonestimable statuses, retained taxa/families, and control reuse. All 253
independent statistics fixtures passed before launch. The checker requires the
producer's exact process identity to finish and its service to report success;
source/script/plan hashes are pinned and rechecked. Its plan and launch record
are `metadata/screened_match_balance_readback_plan_20260928.json` and
`metadata/screened_match_balance_readback_launch_20260928.json`.

Audit resources: one CPU, 24 GiB memory, no swap/GPU/paid resources; estimated
0.1–8 hours and under 0.01 GiB new output. Full output validation remains pending.
These are descriptive balance diagnostics, not calibrated duplication effects.

### Screen sensitivity figures queued behind full validation

The figure stage `scripts/plot_screened_match_balance.py` waits for the exact
independent-audit process to finish successfully and verifies its source receipt
before using the results. It will export PNG, PDF, SVG and a source-value TSV for
S45 and S46, both guides, the alignment-E-value policy and all six screens,
requiring both full and confidence-masked alignments to pass. Three panels per
guide distinguish matched balance, shifts from the original target pool, and
additional shifts from the pre-screen matched pool. Each metric has its own
symmetric color range, fixed across scenarios; cell values and per-screen retained
counts are shown. The TSV preserves feature-specific denominators, including
positive-log-distance exclusions. The scenarios remain illustrative sensitivity
analyses, not designated primary tests.

The stage checks all exported rows against source-derived values and the plotted
arrays against the same matrices. Visual review remains required after generation.
Resources: one CPU, 4 GiB memory, no swap/GPU/charges, 1–10 minutes after the audit,
under 0.1 GiB output. Plan and launch identity are recorded in
`metadata/screened_match_balance_figure_plan_20260928.json` and
`metadata/screened_match_balance_figure_launch_20260928.json`.
