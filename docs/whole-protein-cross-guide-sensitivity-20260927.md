# Exact duplicate pairs across phylogenetic guides

Equal aggregate counts under the MAFFT and profile guides did not establish
identical duplication membership. We now match the exact family and lexical
A/B gene pair across guides, retain guide-specific node identifiers, and keep
all pairs absent from either guide. Within each guide, the matched key is
unique for each threshold. Gene orientation was checked across every row.

The modeled-reference inventories contain 20,178 gene pairs in their union:
14,731 occur in both guides, 2,717 only in MAFFT and 2,730 only in profile. These
are counts within the modeled-reference inventory, not all reconciled events.

| Minimum residues / coverage | All references eligible in both guides | Same structural direction across both guides and all settings | Families | Taxa |
| --- | ---: | ---: | ---: | ---: |
| 30 / 50% | 5,631 | 3,825 | 1,711 | 119 |
| 30 / 70% | 3,436 | 2,358 | 1,108 | 113 |
| 30 / 90% | 930 | 684 | 389 | 81 |
| 50 / 50% | 5,608 | 3,805 | 1,699 | 119 |
| 50 / 70% | 3,434 | 2,356 | 1,106 | 113 |
| 50 / 90% | 930 | 684 | 389 | 81 |

Family and taxon counts describe the same-direction subset. Eligibility requires
three distinct models per reference and passing all 32 mask/order/mapping
configurations. Direction uses the numerical tolerance documented in
[reference sensitivity](whole-protein-reference-sensitivity-20260927.md), not a
biological effect-size cutoff. All six screens are retained; none is selected
because it yields more candidates.

At 50 residues and 70% coverage, 1,411 of the 2,356 pairs also have matching
sequence-divergence and structural-distance directions under both guides.
The other pairs require inspection of the exported sequence-direction states;
they must not collectively be called sequence–structure decoupling. Correlated
signs alone are not a fitted coupling model, and these counts do not account
for dependence among families, taxa or reused structures.

The full outer join retains 121,068 pair/threshold rows. An independent dictionary
join checked every source field, membership status, complete-reference flag and
structural-direction agreement; serialized readback checked every exported cell.
An initial in-memory merge attempt encountered pandas' categorical missing-value
restriction and produced no output; converting the membership indicator to text
before filling missing guide fields resolved it.

Reproduce with `scripts/compare_whole_protein_duplicate_guides.py` (one CPU;
observed runtime under a minute). Large tables and receipt are in
`results/structural_comparisons/whole-protein-cross-guide-sensitivity-20260927-v1`.
Versioned hashes and counts are in
`metadata/whole_protein_cross_guide_sensitivity_completed_20260927.json`.
These results measure reproducibility of exact modeled gene-pair contrasts;
biological duplication effects, phylogenetic uncertainty and independent
validation remain outstanding.
