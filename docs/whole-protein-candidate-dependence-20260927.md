# Dependence in whole-protein duplication candidates

The cross-guide candidate counts do not represent independent biological
observations. We mapped every exact gene pair to all linked A, B and reference
model accession/version identifiers across both guides, retaining all tied
references. For every one of the six coverage screens, we examined both the
complete-reference eligible set and its direction-stable subset.

At the 50-residue/70%-coverage screen:

| Quantity | All-reference eligible | Direction stable |
| --- | ---: | ---: |
| Gene pairs | 3,434 | 2,356 |
| Families | 1,401 | 1,106 |
| Taxa | 122 | 113 |
| Distinct model accession/version combinations | 10,154 | 7,005 |
| Models used by multiple pairs | 153 | 68 |
| Connected groups from shared models | 3,332 | 2,312 |
| Largest group from shared models | 3 | 2 |
| Largest family contribution | 78 | 62 |
| Largest taxon contribution | 401 | 251 |

Adding family membership to the model-sharing graph yields 1,401 and 1,106
connected groups, respectively. In these selected subsets model reuse does not
connect different families. This does not rule out cross-family links elsewhere
in the full inventory or prove statistical independence between families.
Shared ancestry, taxa and prediction errors still need separate treatment.

Both union-find and a full graph traversal independently recovered exactly the
same component memberships for all twelve screen/subset combinations. Complete
serialized readback checked the component table, family/taxon counts and summary.
`component_membership.tsv` retains each gene pair and a deterministic identifier
for its connected group. `sampling_counts.tsv` retains every family and taxon
count rather than only the largest contributors.

Run `scripts/measure_whole_protein_candidate_dependence.py` with one CPU. Outputs
and checksums are in
`results/structural_comparisons/whole-protein-candidate-dependence-20260927-v1`;
the versioned receipt is
`metadata/whole_protein_candidate_dependence_completed_20260927.json`.
These are dependencies for the statistical analysis, not corrected significance
tests or an estimate of an effective sample size. See
[cross-guide selection](whole-protein-cross-guide-sensitivity-20260927.md).
