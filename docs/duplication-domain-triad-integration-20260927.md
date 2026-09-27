# Domain comparisons integrated with duplication events and sequence divergence

The normalized SQLite dataset links the existing provisional duplicate/reference
triads to domain comparison coverage and sequence-tree covariates. It retains
all reference ties and both species-tree guides. These are dependent alternatives,
not extra independent observations. A provisional reference is not automatically
a biological ortholog or ancestral structure.

Inputs are frozen in
[the integration plan](../metadata/duplication_domain_triad_integration_plan_20260927.json).
The source tables contain 147,776 triad/policy rows (36,944 event/reference
combinations under four annotation policies), 301,380 domain policy/boundary
links, 140,790 domain pair/mask coverage rows, and 34,909 guide/event sequence
covariate rows. Data remain outside Git history.

The database contains these tables and views:

| Relation | Purpose |
|---|---|
| `triads` | All event/reference/policy records, including identical-model and incomplete-annotation cases |
| `domain_links` | Protein pair and Pfam/policy/boundary to exact domain intervals and domain-pair IDs |
| `coverage` | Both-order metrics, six coverage screens and exclusions for every domain pair/mask |
| `sequence_covariates` | Original gene-oriented sequence-tree measurements and chosen reference |
| `matched_domains` | Pfam candidates present in all three comparisons, with identical intervals at each shared protein |
| `domain_availability` | Every triad/policy under both boundaries, including explicit absence of a common candidate |
| `triad_domain_coverage` | Each matched domain under both masks, six all-three-comparisons screen flags and sequence covariates |
| `provenance` | Exact source and script hashes |

The three comparisons are duplicate A–B, A–reference and B–reference. All
must match the same Pfam candidate under the same policy and interval-boundary
choice. The A interval in A–B must equal that in A–reference; the same rule
applies to B and the reference. Any inconsistent shared vertex stops the build.
The producer checks that no matched row is lost during coverage/sequence joins.

A coverage screen passes a triad/domain/mask only if all three comparisons
pass; each comparison already requires both input orders. No favorable order,
reference, guide, boundary or policy is chosen, and no result is averaged.
The numerically discrepant pair/mask remains quarantined by the source screen.

Gene-oriented tip divergence and duplicate-pair distance retain their source
reference identifier. Distances to **this particular reference** are populated
only when its gene ID exactly matches the reference used for the sequence
calculation; otherwise they are SQL NULL with `different_tied_reference`.
A tied alternative therefore never inherits another protein's reference-specific
distance. These are whole-protein sequence-tree covariates, not domain-specific
sequence rates or calendar-time rates.

The independent checker compares all imported CSV fields, reconstructs every
three-way Pfam intersection and shared-vertex interval match in Python, checks
the complete availability and matched sets, then checks every coverage flag,
sequence field and reference status against source rows. It also reports unique
guide/event keys so policy and domain multiplicity cannot be mistaken for
independent duplication events.

Reproduce with new output paths in a new plan when previous artifacts exist:

```bash
python scripts/integrate_duplication_domain_triads.py --plan metadata/duplication_domain_triad_integration_plan_20260927.json
python scripts/check_duplication_domain_triad_integration.py --plan metadata/duplication_domain_triad_integration_plan_20260927.json --output metadata/duplication_domain_triad_integration_completed_20260927.json
```

This integration does not replace the failed strict alignment audit, establish
coordinate rank, validate domain boundaries or infer a duplication effect.
`biological_inference_eligible` remains false. Missing common candidates are
annotation/comparison availability outcomes, not evidence of domain loss.
Subsequent analyses need explicit geometry checks, appropriate control groups,
phylogenetic and family dependence, reference/policy sensitivity, and selection
of scientifically interpretable structural summaries.

## Verified integration results

The full independent readback passed: 72,336 matched domain-triad records,
144,672 records after expanding both masks, and 295,552 complete
triad/policy/boundary availability records. The matched records cover 6,851
guide/event keys across the two alternative guides; this total must not be
interpreted as 6,851 independent biological events.

For the 30-residue/70%-coverage screen under pLDDT70, at least one complete
three-comparison domain record passes for 3,127 MAFFT-guide events and 3,129
profile-guide events. These are availability counts across any retained
reference/policy/boundary alternative, not counts stable across all alternatives,
not additional observations to pool, and not evidence of structural asymmetry.

Of the 144,672 domain/mask records, 134,720 use the exact sequence-covariate
reference and 9,952 use another tied reference. The latter retain explicit NULLs
for distances to that particular reference. The
[completion record](../metadata/duplication_domain_triad_integration_completed_20260927.json)
provides all six screens and exact source/readback hashes.
