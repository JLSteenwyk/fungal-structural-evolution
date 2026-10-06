# Provisional hydrophobin profile inventory

This is a reproducible candidate inventory for a later, carefully controlled
structural case study. It is not an orthology call, a hydrophobin functional
assignment, a structural clade, a duplication count, or an evolutionary result.

The inventory scans all 64 completed shards of the full-Pfam annotation catalog
for exact matches to the three current Pfam hydrophobin profiles. The catalog
metadata records `PF01185.24` (Hydrophobin) and `PF06766.17` (Hydrophobin_2) as
Families, and `PF29785.1` (Hydrophobin_D) as a Domain. The latter distinction is
preserved rather than converting every profile into a family call.

The immutable receipt is
[`metadata/hydrophobin_profile_inventory_20261006_v1.json`](../metadata/hydrophobin_profile_inventory_20261006_v1.json).
It records 1,292 accepted profile hits, yielding 1,318 taxon--protein--profile
rows across 1,290 distinct sequences and 220 of the 526 analysis-manifest taxa.
There are 1,064 rows for `PF01185.24`, 253 for `PF06766.17`, and one for
`PF29785.1`. The row total is greater than the hit total because a sequence can
have more than one retained protein/taxon linkage in the full domain database.

Of the 1,318 rows, 803 link to a record in the current source-structure
registry and 515 do not. This is source-model availability only: it has not yet
checked coordinate retrieval, residue confidence, domain completeness, PAE,
or suitability for a pairwise structural analysis.

An independent reconstruction in
[`metadata/hydrophobin_profile_inventory_readback_20261006_v1.json`](../metadata/hydrophobin_profile_inventory_readback_20261006_v1.json)
rescans all three profile-bearing catalog shards, rebuilds every protein/taxon
link from the full-domain database and rechecks source-model availability in
the independently read-back structure registry. It reproduces all 1,292 hits
and 1,318 rows exactly. This establishes inventory provenance only; it does
not upgrade candidates to homologs or structural observations.

## Reproduction

The result directory remains outside Git because it contains the full candidate
table. The receipt records its SHA-256 digest. The command rejects an existing
output or published receipt so that a result cannot be silently overwritten.

```bash
python scripts/inventory_hydrophobin_candidates_v1.py \
  --catalog results/domains/full-annotation-catalog-v1/annotation_shards.tsv \
  --domain-db results/domains/full-domain-database-v1/domains.sqlite \
  --domain-receipt metadata/full_domain_database_receipt.json \
  --structure-registry results/domains/whole-proteome-structure-domain-registry-20261005-v1/structure_domains.sqlite \
  --registry-receipt results/domains/whole-proteome-structure-domain-registry-20261005-v1/receipt.json \
  --manifest metadata/analysis_manifest.tsv \
  --output results/domains/hydrophobin-profile-inventory-20261006-v1 \
  --published-receipt metadata/hydrophobin_profile_inventory_20261006_v1.json

python scripts/readback_hydrophobin_candidates_v1.py \
  --catalog results/domains/full-annotation-catalog-v1/annotation_shards.tsv \
  --domain-db results/domains/full-domain-database-v1/domains.sqlite \
  --structure-registry results/domains/whole-proteome-structure-domain-registry-20261005-v1/structure_domains.sqlite \
  --manifest metadata/analysis_manifest.tsv \
  --inventory results/domains/hydrophobin-profile-inventory-20261006-v1/hydrophobin_candidates.tsv \
  --producer-receipt metadata/hydrophobin_profile_inventory_20261006_v1.json \
  --receipt metadata/hydrophobin_profile_inventory_readback_20261006_v1.json
```

The next gate is a family-specific analysis that first establishes homologous
groups and reconciled gene trees, then qualifies coordinates and uncertainty,
compares domains with direct geometry and structural alphabets, controls for
prediction source, and only then assesses sequence--structure change. Any
biological interpretation also requires appropriate functional evidence and
validation.
