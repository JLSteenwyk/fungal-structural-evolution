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
