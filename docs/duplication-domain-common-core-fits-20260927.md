# Structural distances on identical residue triples

This stage fits all 7,986 oriented domain triads under both masks, all eight
native alignment-order combinations and both common-residue definitions:
255,552 output rows. Identical mappings may reuse a cached numerical fit;
all output alternatives and links to the 72,336 source event/domain records
remain explicit. No input order or mapping definition is selected for favorable
results.

For each core with at least three residues and no upstream edge exclusion,
three independent least-squares proper rigid fits give RMSDs for A–B,
A–reference and B–reference. Every fit uses the **same residue triples** and
therefore the same count. Refitting does not impose a single global rotation
on all three proteins; each pair minimizes its own coordinate discrepancy.
The signed quantity `rmsd_ar_minus_br` preserves duplicate A/B orientation.
It is a difference in reference similarity, not a branch length, ancestral
change, calendar-time rate, or proof of asymmetric evolution.

The two core definitions are:

- `reference_common`: residues shared by the two reference-based mappings.
- `cycle_consistent`: the subset also agreeing with the direct duplicate-pair
  correspondence.

A lack of consistency around the three maps is preserved as
`mapping_disagreement_count`; it does not automatically discard a nonempty
reference-common core. Numerical discrepancies, unavailable inputs and
original degenerate-edge fits remain source exclusions for both definitions.
Such rows are retained with blank fitted values. Cores shorter than three
residues also retain blank fits. Longer cores with degenerate refitted geometry
retain descriptive values and an explicit failure flag.

All three core fits are checked for proper-rotation uniqueness using the
previously verified geometry method. The six length/coverage screens are
reapplied **after** intersection: at least 30/50 common residues and 50/70/90%
coverage of each original interval. Coverage denominators never shrink with
confidence masking or intersection. Screen flags require all three core fits
to be numerically unique and all source exclusions to be absent. These are
exploratory technical screens, not validated biological thresholds.

The table also reports exact pairwise sequence identity on the fitted triples,
mean PDB pLDDT for each protein, and the fraction of triples where all three
residues have pLDDT ≥70. These confidence summaries use the serialized PDB
values and do not calibrate structural error or validate alignment homology.

Fixtures check same-core RMSD metric bounds, duplicate-label reversal, rigid
coordinate transformations, exact sequence identities, confidence summaries,
collinear degeneracy and rejection of two-residue fits. The full producer
still requires an independent numeric readback before downstream use.

The [pinned plan](../metadata/duplication_domain_common_core_fits_plan_20260927.json)
allows one CPU, 16 GiB memory, no swap and 2 GiB output, with an initial 2–30
minute planning range. No GPU prediction, native alignment or external charges
are involved. Outputs refuse overwrite.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/check_common_residue_domain_fit_cases.py
OPENBLAS_NUM_THREADS=1 python scripts/fit_domain_triad_common_residues.py --plan metadata/duplication_domain_common_core_fits_plan_20260927.json
```

Phylogenetic/family dependence, reference choice, order/mask robustness,
confidence and alignment uncertainty still need to be carried into subsequent
analyses. A signed distance contrast alone does not establish a duplication
effect or identify which copy changed after duplication.
