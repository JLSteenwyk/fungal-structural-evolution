# Qualified AlphaFold ecological coverage

The expanded AlphaFold paired-input collection was assessed across all 32
species in the frozen sample-linked ecological evidence table. It contains
125 eligible markers overall; 23 curated species have at least one eligible
marker. Of all 496 unordered taxon pairs, 246 share at least one marker with
50 or more jointly observed AA/3Di alignment columns. These are descriptive
coverage counts, not independent ecological contrasts or a power calculation.

The existing confidence masks and eligibility rules are preserved. All taxon,
pair, pair-marker and provisional-group counts passed independent verification
using separate FASTA parsing and integer mask matrix multiplication. The source
paired inputs also have a completed full reconstruction from qualified arrays.

| Provisional group | Taxa | AlphaFold markers with ≥50 columns in every member | ESMFold counterpart |
|---|---:|---:|---:|
| Amanita | 6 | 0 | 31 |
| Cantharellales | 4 | 0 | 0 |
| Cenococcum and comparators | 3 | 1 | 60 |
| Laccaria | 2 | 99 | 59 |
| Phallomycetidae | 4 | 0 | 35 |
| Piloderma | 1 | 119 | 63 |
| Suillus | 9 | 0 | 47 |
| Thelephora | 2 | 0 | 0 |
| Tuber | 1 | 115 | 62 |

ESMFold counts come from the separately verified
[September 22 dataset](qualified-ecology-coverage-20260922.md). The table compares
source-specific coverage; it does not merge observations, assume identical
retained columns across sources or measure agreement between predictions.
Singleton and same-state groups do not supply ecological transitions.

AlphaFold provides 112 qualified shared markers (31,664 common columns) for
Hydnum rufescens–Botryobasidium botryosum, and 96 (25,321 columns) for
Hydnum–Tulasnella calospora. Tulasnella remains an orchid-mycorrhizal species,
not an absence-of-symbiosis control. Cantharellus anzutake lacks eligible
AlphaFold markers; Hydnum, Botryobasidium and Tulasnella lack eligible markers
in the frozen ESMFold dataset. Neither source therefore provides the complete
four-taxon Cantharellales group.

The two Thelephora species likewise lack complete common-source coverage.
Amanita and Cenococcum/comparator coverage is stronger in ESMFold. These results
identify which source can support further comparison, but do not resolve
orthology, ecological-state uncertainty, ancestral direction or independent
transition replication. The prior phylogenetic classification diagnostic does
not establish independent origins; ecological association testing remains open.

## Reproduction and evidence

```bash
python scripts/assess_qualified_afdb_ecology_overlap.py --plan metadata/qualified_afdb_ecology_overlap_plan_20260927.json
python scripts/readback_qualified_ecology_overlap.py --plan metadata/qualified_afdb_ecology_overlap_plan_20260927.json --output metadata/qualified_afdb_ecology_overlap_readback_20260927.json
```

Use a copied plan and fresh output paths for reruns. The AlphaFold producer is
a separate variant, preserving the pinned ESMFold implementation and adding an
explicit assertion of the 50-column threshold used in output field names.
The original independent matrix checker works without modification. The plan
pins inputs, proofs and scripts. Full tables are outside Git under
`results/ecology/qualified-afdb-recovered-overlap-20260927-v1`; compact taxon and
group tables, receipt and independent proof are archived as
`metadata/qualified_afdb_ecology_*_20260927.*`.

The pre-run allowance was one CPU, 4 GiB memory, 0.02 GiB output and 1–10 minutes;
production and readback completed in seconds. No GPU predictions or paid
resources were used. The full project's ecological aim remains unfinished.
