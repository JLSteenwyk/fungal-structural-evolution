# Sensitivity of whole-protein duplicate/reference contrasts

The common-core fits provide signed structural contrasts
`RMSD(A, reference) − RMSD(B, reference)` and sequence-divergence contrasts
`identity(B, reference) − identity(A, reference)`. Positive values mean A is
farther from the reference in the respective measurement. A/B labels retain
the recorded gene ordering; the reference is an extant protein, not an ancestor.

`scripts/summarize_whole_protein_reference_sensitivity_v2.py` checks source hashes
and reconstructs all triad contrast ranges from 563,808 fit rows. It uses both
confidence masks, all eight alignment-order combinations and both common-residue
mapping definitions. A second full CSV pass independently checks every range.

For each of six length/coverage screens, a reference is eligible only if it has
three distinct models and passes all 32 configurations. Event ranges encompass
all eligible references. Every event is retained, including those without an
eligible reference. `all_references_eligible` distinguishes complete reference
coverage from cases where only a subset qualifies; a stable sign in the latter
case does not establish robustness to the missing references. Every reference
link is also exported with its unfiltered range and coverage flags.

The numerical-zero tolerance is 1e-8, above the maximum observed independent
fit discrepancy. It is not a meaningful biological effect-size threshold.
Ranges strictly above or below this tolerance receive directional labels;
values entirely within the tolerance are numerically zero; other ranges are
setting/reference sensitive. Six boundary cases verify this classification.
This step does not fit a statistical model or generate p-values. Variation
across analysis choices is not prediction uncertainty or a confidence interval.

These descriptive results prepare reference-robustness and sequence/structure
comparisons for subsequent phylogenetic analysis. Family dependence, reused
models, shared taxa, guide-tree differences, domain orientation, controls and
prediction-source uncertainty still require assessment before interpreting
asymmetric evolution or excess structural change.

The first export attempt stopped during readback because pandas parsed numeric
columns containing blanks as strings. The v1 partial output is preserved and
has no completion receipt. Version 2 checks blank locations exactly, compares
finite range values numerically, and compares other fields as strings. This
changes serialization validation only; input selection and range calculations
are unchanged.

The corrected full run completed 17,619 triad ranges, all 36,944 event/reference
links and 209,454 event/screen rows. Source-range reconstruction and serialized
readback passed. Outputs are in
`results/structural_comparisons/whole-protein-reference-sensitivity-20260927-v2`;
provenance and hashes are in
`metadata/whole_protein_reference_sensitivity_completed_20260927.json`.

At the 50-residue/70%-coverage screen, each guide has 3,534 events with at least
one eligible reference. Of these, 3,530 have **all** tied references eligible:
1,292 have a consistently positive structural contrast, 1,131 consistently
negative, and 1,107 a sign sensitive to the settings or references. The remaining
four events have only partial reference eligibility. These numerical sign
counts do not establish a biologically meaningful magnitude, independent
replication or significant asymmetry. Identical guide-level counts do not prove
identical event membership.
