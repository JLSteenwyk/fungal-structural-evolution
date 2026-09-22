# Ecological sampling in the confidence-filtered ESMFold alignments

The completed 25,322-model ESMFold integration produces eligible paired
AA/3Di alignments for 122 markers. Of the 32 species with curated ecological
classifications, 28 have at least one eligible marker. All 496 species pairs
were evaluated; 378 have at least one marker with 50 or more common qualified
alignment columns. This is a coverage diagnostic, not 378 independent
contrasts or evidence of an ecological effect.

This analysis intersects the actual observed columns of both species after
the existing confidence and taxon/marker eligibility rules. Sharing a model
or a marker name is insufficient. The AA and 3Di masks must agree exactly.
The descriptive 50-column overlap threshold does not establish statistical
power and is not an additional fitted-analysis eligibility rule.

| Provisional comparison group | Taxa | Markers with >=50 columns observed in every member |
|---|---:|---:|
| Amanita | 6 | 31 |
| Cantharellales | 4 | 0 |
| Cenococcum and comparison species | 3 | 60 |
| Laccaria | 2 | 59 |
| Phallomycetidae | 4 | 35 |
| Piloderma | 1 | 63 |
| Suillus | 9 | 47 |
| Thelephora | 2 | 0 |
| Tuber | 1 | 62 |

Singleton groups describe coverage only. Groups were inherited from the
curation table, not inferred as independent transitions here. Same-state
groups such as Laccaria and Suillus do not themselves provide ecological
contrasts. The Phallomycetidae group retains the source-specific Ramaria
rubella saprotrophic classification and its previously documented naming
and genus-classification caveats.

Hydnum rufescens, Botryobasidium botryosum, Tulasnella calospora and Thelephora
terrestris have no eligible markers in this ESMFold-only dataset. This is not
absence of public models or proof of biological absence. Refreshed AlphaFold
processing remains separate; combining prediction sources requires explicit
source controls, and does not automatically solve common-source coverage.

Reproduce the full analysis using a fresh output directory in a copied plan:

```bash
python scripts/assess_qualified_ecology_overlap.py \
  --plan metadata/qualified_esmfold_ecology_overlap_plan.json
python scripts/readback_qualified_ecology_overlap.py \
  --plan metadata/qualified_esmfold_ecology_overlap_plan.json \
  --output results/ecology/qualified-overlap-readback-new.json
```

The pinned plan checks the ecological evidence receipt, paired-input receipt,
independent array readback and every alignment artifact. The full output is
`results/ecology/qualified-esmfold-overlap-20260922-v1`, including all pair and
pair-marker counts. Compact taxon/group tables and receipts are archived in
`metadata/qualified_esmfold_ecology_*`. Independent manual FASTA parsing and
integer mask matrix multiplication reproduce every count across all taxa,
markers, pairs and groups. This does not repeat confidence calibration or
establish biological homology.

The pre-run resource allowance was one CPU, 4 GiB memory, 0.02 GiB output and
1–10 minutes; the analysis and independent check each completed in seconds.
No GPU prediction or paid infrastructure was used. Phylogenetic transition
review, predictor sensitivity and ecological effect testing remain required.
