# Cross-guide identities and sequence/structure direction

This comparison uses an exact outer join on protein family, lexical gene A,
gene B, Pfam accession and coverage screen. Guide-specific gene-tree node names
are retained but are not treated as cross-guide identifiers. A repeated event
for the same exact join key within a guide stops the producer. Shared pairs
must have the same taxon identity. All six screens and three descriptive
margins remain in the output.

Structural direction refers to the range of common-core RMSD(A, reference) minus
RMSD(B, reference), already requiring every retained reference, annotation
policy/boundary and fitting alternative. Sequence direction refers to the
previously checked whole-protein sequence-tree tip divergence of the same
lexically oriented gene pair. A sequence direction is usable only when the
pair-distance status is resolved and the direction is A-longer or B-longer.
Unresolved, incomplete, missing-guide and guide-sensitive relationships remain
explicit. No signed distance is relabeled as an ancestral change.

## Independently verified results

At 30 common residues, 70% original interval coverage and the descriptive
0.1 Å margin, the union has 4,320 event/domain identities: 4,143 shared by both
guides, 90 profile-only and 87 MAFFT-only. Of the shared identities, 947 retain
the same stable structural direction in both guides and 3,196 are unresolved
in one or both. None has opposite stable structural directions at this screen
and margin. The 947 combinations represent 894 distinct family/gene-pair keys;
multiple domains within a pair must not be counted as independent events.

Among those 947 shared stable structural combinations:

| Relationship to whole-protein sequence direction | Event/domain combinations |
|---|---:|
| Concordant in both guides | 739 |
| Discordant in both guides | 202 |
| Sequence direction unresolved | 6 |

These are descriptive concordance counts, not a phylogenetically adjusted
association test. Both guide analyses reuse the same biological observations
and predictions. The structures derive from sequence, the reference is not an
ancestor, and sequence covariates are whole-protein rather than domain-specific.
The 202 discordant combinations are candidates for closer examination; they do
not establish structural acceleration, decoupling, adaptive evolution or a
functional difference.

The independent checker reconstructed the complete outer join through dataframe
merges, compared every source field and key, and independently calculated every
directional label/count across all 77,760 comparison rows from 50,778 source
rows. Missing guide values remain blank rather than becoming zeros. See the
[completion record](../metadata/duplication_domain_guide_comparison_completed_20260927.json)
and [pinned plan](../metadata/duplication_domain_guide_comparison_plan_20260927.json).

```bash
python scripts/compare_domain_robustness_guides.py --plan metadata/duplication_domain_guide_comparison_plan_20260927.json
python scripts/check_domain_guide_comparison.py --plan metadata/duplication_domain_guide_comparison_plan_20260927.json --output metadata/duplication_domain_guide_comparison_completed_20260927.json
```

Existing outputs refuse overwrite. This stage performs table comparisons only,
with no GPU, new structure inference or external charges. Phylogenetic/family
controls, taxonomic sampling assessment, matched controls, prediction-source
sensitivity and functional case review remain necessary.
