# Common residue correspondences for duplicate/reference domains

Pairwise RMSDs may use different residue sets, so their difference can confound
structural divergence with alignment coverage. This stage establishes explicit
common residue sets before calculating any structural contrast.

The fully checked domain-triad integration contains 72,336 event/domain records
but only 7,986 distinct oriented interval triples (duplicate A, duplicate B,
reference). Each is evaluated under full and pLDDT70 masks and all eight
combinations of native input orders for A–B, A–reference and B–reference:
127,776 dispositions. Computational reuse preserves every guide, family, gene
node, gene pair, tied reference, annotation policy, boundary and Pfam link.
These links remain dependent alternatives, not independent observations.

For each disposition:

1. Read the three hashed native checkpoints and convert their aligned residue
   indices to original protein residue positions using the frozen input manifest.
2. Normalize mapping direction to A→B, reference→A and reference→B, irrespective
   of native input order.
3. Intersect reference positions in the two reference-based mappings. Retain
   the corresponding triples `(A_position, B_position, reference_position)`.
4. Separately retain the subset whose A→B correspondence agrees with the direct
   duplicate-pair alignment. This reports consistency around all three mappings;
   it does not prove that they are true evolutionary homologies.
5. Preserve every original numerical/geometry/input exclusion and explicitly
   mark missing common-reference or consistent triples. Mapping availability
   never clears a pre-existing quarantine.

The output separates `reference_common_triples` from
`cycle_consistent_triples`; it does not hide disagreement by relabeling the
intersection. Original interval lengths are carried for subsequent coverage
checks. Counts are occurrences across masks/order combinations and must not
be treated as unique residues, proteins, events or independent measurements.
There is no order averaging or selection. Later work must evaluate robustness
across orders and masks, recompute fits on the same residue triples, and account
for loss of coverage when forming intersections.

This stage performs no new structural alignment, coordinate fitting, structure
prediction or inference of biological asymmetry. It is conditional on the
frozen domain-boundary and alignment hypotheses. Independent verification of
the complete mapping inventory remains required before downstream fitting.

The [plan](../metadata/duplication_domain_common_residues_plan_20260927.json)
pins the completed geometry readback and domain-triad integration. It allows
one CPU, 16 GiB RAM, no swap, and 4 GiB output; the initial runtime planning
range is 1–20 minutes. No external charges or GPU use are involved.

```bash
python scripts/check_domain_common_residue_cases.py
python scripts/prepare_domain_triad_common_residues.py --plan metadata/duplication_domain_common_residues_plan_20260927.json
```

Outputs refuse overwrite. Fixtures check common-reference intersections,
inconsistent direct mappings, disjoint and missing maps, noncontiguous original
positions, nonbijective-map rejection and swapping duplicate labels.
