# Structural coverage of the 45-species ecological evidence set

Updated both prediction-source coverage analyses using the complete reviewed
45-species evidence table. Each analysis evaluated all 990 unordered pairs,
every available marker, and the inherited nine provisional groups plus thirteen
unassigned singleton labels. Independent FASTA parsing and matrix multiplication
verified every taxon, pair, pair-marker and group count. The singleton labels
are bookkeeping identifiers, not established ecological transitions.

| Coverage measure | AlphaFold | ESMFold |
|---|---:|---:|
| Eligible markers in the source collection | 125 | 122 |
| Curated taxa with any eligible marker | 34 | 32 |
| Pairs sharing a marker with ≥50 qualified columns | 545 | 496 |

All 45 taxa have an eligible marker in at least one source: 21 in both, 13 only
in AlphaFold and 11 only in ESMFold. For pair coverage, 195 pairs qualify in
both source datasets, 350 only in AlphaFold, 301 only in ESMFold and 144 in
neither. Thus 846 pairs meet the descriptive shared-column criterion in at
least one source. These are coverage counts, not independent contrasts,
statistical power estimates or an ecological association result.

“Both sources” does not assert identical marker membership or homologous
retained columns across datasets. All intersections were calculated within a
source's actual paired AA/3Di masks. No cross-source structural observations
were pooled. More complete acquisition does not remove source-dependent
missingness, uncertainty in species labels, or shared ancestry.

Every original row from the 32-species analyses was recovered exactly in each
expanded result: 32 taxon rows, 496 pair rows and nine group rows. Their reported
Amanita, Cantharellales, Cenococcum and other group limitations therefore remain
unchanged. The new source classifications are not automatically negative ECM
observations; ecological trait modeling still requires explicit state coding
and independently supported transitions.

## Reproduction

```bash
python scripts/assess_qualified_afdb_ecology_overlap.py --plan metadata/qualified_afdb_reviewed45_ecology_plan_20260927.json
python scripts/assess_qualified_ecology_overlap.py --plan metadata/qualified_esmfold_reviewed45_ecology_plan_20260927.json
python scripts/readback_qualified_ecology_overlap.py --plan metadata/qualified_afdb_reviewed45_ecology_plan_20260927.json --output metadata/qualified_afdb_reviewed45_ecology_readback_20260927.json
python scripts/readback_qualified_ecology_overlap.py --plan metadata/qualified_esmfold_reviewed45_ecology_plan_20260927.json --output metadata/qualified_esmfold_reviewed45_ecology_readback_20260927.json
python scripts/compare_reviewed_ecology_source_coverage.py
```

Use fresh output paths for reruns. Full tables remain outside Git under
`results/ecology/qualified-{afdb,esmfold}-reviewed45-overlap-20260927-v1`.
Source-bound plans, receipts and full independent proofs are archived under
`metadata/qualified_{afdb,esmfold}_reviewed45_ecology_*_20260927.json`.
The comparison script checks artifact hashes, matches species and state fields,
checks preservation of old outputs and records source-specific coverage in
`metadata/reviewed45_ecology_source_coverage_20260927.tsv` and its receipt.

Each source analysis retained its one-CPU, 4-GiB, 1–10-minute planning allowance
and completed in seconds; no GPU prediction or paid infrastructure was used.
This advances the ecological coverage requirement without completing aim 5.
