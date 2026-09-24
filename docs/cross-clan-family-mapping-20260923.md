# Cross-clan candidate family mapping

The full inventory of 68 boundary/cluster/Pfam-pair candidates now has source
protein, taxon, and family assignments under both the profile and MAFFT guides.
The 4,600 candidate-member rows expand to 4,746 protein-linked rows, representing
2,759 distinct source proteins across 224 taxa. All members are retained, before
confidence or direct-alignment screening; boundary views overlap.

Results are outside Git in
`results/structural_clusters/cross-clan-family-mapping-20260923-v1/`:

- `candidate_proteins.tsv`: complete member-to-protein identities and both family assignments.
- `candidate_family_summary.tsv`: 136 rows covering all 68 candidates, separately
  for all members and pair/cluster-exclusive models. Counts include each side
  and the intersection for models, proteins, taxa, and both family assignments.
- `receipt.json`: source hashes, output hashes, scope, and counts; archived as
  `metadata/cross_clan_family_mapping_completed_receipt_20260923.json`.

Among exclusive-model summaries, two candidate entries share a family between
the two Pfam sides under each guide. This is descriptive assignment overlap,
not proof of homology, orthology, duplication, or a structural transition.
Subsequent case selection must join the direct-comparison confidence results,
inspect domain architectures, and evaluate the corresponding gene trees.

Validation independently rebuilt the complete export using a SQL join and
compared all 4,746 serialized identities. Pandas/set aggregation checked every
count in all 136 summary rows. An initial audit query was stopped because its
query plan scanned the candidate table repeatedly; explicit join order completed
the same full validation. No analysis producer was interrupted.

Reproduce with one CPU (observed runtime seconds, excluding the inefficient
initial audit), no GPU, and local input databases:

```bash
python scripts/map_cross_clan_families.py \
  --candidates results/structural_clusters/cross-clan-candidate-inventory-20260923-v1 \
  --composition results/structural_clusters/domain-cluster-composition-20260923-v1 \
  --composition-audit metadata/domain_cluster_composition_completed_readback_20260923.json \
  --output results/structural_clusters/cross-clan-family-mapping-20260923-v1
```

The output directory must not already exist.

## Direct-comparison join within shared families

`scripts/screen_cross_clan_shared_families.py` joins every exclusive-model
cross-clan pair within an assigned family to the completed bidirectional
structural comparisons. The two qualifying candidate entries are the alignment
and envelope views of PF00085.27 versus PF26973.1 in the same structural cluster.
Both guides assign the shared family as OG0000030.

The complete join contains 1,312 protein-pair/guide/boundary rows representing
631 distinct interval pairs. These repeated views and source-protein links are
not independent observations. For each guide separately, rows passing minimum
TM-score 0.5 and minimum coverage 0.8 are:

| Boundary | No confidence fraction cutoff | 0.5 | 0.8 | 0.9 |
|---|---:|---:|---:|---:|
| Alignment | 172 | 172 | 172 | 169 |
| Envelope | 375 | 371 | 367 | 342 |

Confidence fractions require pLDDT at least 70 in both whole domains and matched
residues in both alignment directions. No PAE qualification is included. These
screens prioritize gene-tree and domain-architecture inspection; they are not
tests of evolutionary transitions or homology.

All serialized pair identities were independently reconstructed with pandas
merges, all direct scores matched their source table, and every screen flag was
recomputed. Results reside in
`results/structural_comparisons/cross-clan-shared-family-screen-20260923-v1/`;
the source/output hash receipt is archived in
`metadata/cross_clan_shared_family_screen_completed_receipt_20260923.json`.
The join uses one CPU, takes seconds, and requires no new structural inference.

```bash
python scripts/screen_cross_clan_shared_families.py \
  --mapping results/structural_clusters/cross-clan-family-mapping-20260923-v1 \
  --comparisons results/structural_comparisons/cross-clan-confidence-summary-20260923-v1 \
  --output results/structural_comparisons/cross-clan-shared-family-screen-20260923-v1
```
