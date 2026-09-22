# Complete local prediction inventory checkpoint

The five completed ESMFold cohorts contain 25,322 distinct exact protein
sequences and model IDs. The combined inventory preserves every original
model record unchanged, including coordinates, prediction receipt, PAE file
binding and prediction configuration. It includes the complete 10,522-model
original cohort and excludes its overlapping earlier partial snapshot.

| Cohort | Models | Configurations |
|---|---:|---:|
| Original short-protein queue | 10,522 | 1 |
| Follow-on short-protein queue | 4,252 | 1 |
| Ecology-targeted short-protein queue | 675 | 1 |
| 513–768 residues | 5,510 | 1 |
| 769–1,024 residues | 4,363 | 3 |

The seven recorded configurations differ only in input-receipt hash, maximum
sequence length, and physical GPU identity. All complete configurations and
per-model batch assignments remain explicit. This observation does not prove
the absence of batch effects or establish comparable biological sampling.

`scripts/combine_completed_prediction_inventories.py` built the inventory
under `results/structures/esmfold-all-completed-20260922-v1/` in 20 seconds.
It verified source artifact hashes, every coordinate and prediction receipt,
per-model configuration binding, exact counts and disjoint sequence/model IDs.
It preserves PAE hashes from the original records; it does not repeat the
independent PAE or coordinate-accuracy audits. The pre-launch plan was one
CPU equivalent, 8 GiB RAM, no swap, 2 GiB output, a 50 GiB free-disk gate and
0.1–4 hours of uncalibrated runtime planning on existing local resources.

`scripts/readback_completed_prediction_inventory.py` independently compared
every combined record, cohort assignment and complete configuration with its
sources. The full readback passed. Receipts are archived as
`metadata/esmfold_all_completed_inventory_receipt.json` and
`metadata/esmfold_all_completed_inventory_readback.json`.

Exact-sequence matching to the global marker inventory yields 25,509
taxon/marker/protein links across 295 taxa and 124 markers. The reference
inventory contains 58,883 unique sequences and 59,840 links. These numbers
describe this local ESMFold inventory alone, before confidence qualification;
they do not include separate AlphaFold coverage and are not whole-proteome
coverage estimates. Taxon, marker and cohort tables and the full matching-link
table are in `results/structures/esmfold-all-completed-readback-20260922-v1/`.

Reproduce the full readback into a fresh output directory with:

```
python scripts/readback_completed_prediction_inventory.py \
  --plan metadata/esmfold_all_completed_inventory_plan.json \
  --global-links results/phylogeny/markers-full-v1/protein_mapping.tsv \
  --output results/structures/esmfold-all-completed-readback-new
```

The full original cohort is undergoing residue mapping; long and extended
cohorts have completed residue mapping and are undergoing feature/confidence
qualification. Combined qualified mappings, regenerated paired alignments,
updated phylogenetic fits and evolutionary tests remain unfinished.
