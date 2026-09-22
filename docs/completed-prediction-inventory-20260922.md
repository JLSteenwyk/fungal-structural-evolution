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

## Combined residue mapping

The combined mapping started under
`fungal-all-completed-mapping-20260922.service`, using
`scripts/map_complete_prediction_union.py` and
`metadata/esmfold_all_completed_mapping_plan.json`. The plan pins 413 source
and script files. It requires the passed inventory readback, binds every
model to its previously audited prediction record, and preserves all original
audit statuses and table schemas in a derived union. A missing source column
is serialized as blank, not inferred from another cohort. A focused check
verified this handling for the ecology audit's absent reuse-state column.

The controller maps all 25,322 models, then independently verifies every
exported residue against the original NPZ sequence/confidence, Stockholm
alignment and frozen profile matrix. The declared scope is all 25,509 exact
global marker/protein links. The mapped snapshot is
`results/structural_markers/esmfold-all-completed-20260922-v1`; its independent
readback is `results/structural_markers/esmfold-all-completed-residue-readback-20260922-v1`.

The pre-launch allowance is one CPU equivalent, 32 GiB RAM, no swap, 50 GiB
output, a 100 GiB free-disk gate, and 1–24 hours of uncalibrated planning. The
larger memory allowance supports full residue-level dictionaries during
readback. This uses existing predictions and does not launch GPU inference.
Completed mapping and confidence-qualified encoding integration remain
separate gates before rebuilding paired phylogenetic inputs.

## Qualified encoding integration

`scripts/merge_complete_qualified_encodings.py` is queued under
`fungal-all-completed-encodings-20260922.service`. Its plan,
`metadata/esmfold_all_completed_encoding_union_plan.json`, waits for four
verified process identities: combined residue mapping and the original, long,
and extended cohort feature controllers. Each predecessor must finish with
the expected successful receipt bound to its pinned plan. Existing follow-on
and ecology qualification receipts and full original-PAE-context readbacks
are pinned directly.

The integration requires all five cohort qualified summaries to match their
mapping and directional PAE readbacks. It checks every model's sequence,
provenance, encoding hash and summary totals, then compares the complete
combined marker-link and residue-link tables against the disjoint multiset
union of the five source mappings. Raw inventory records acquire only the
standard mapping `source_record_id` field; prediction batch provenance stays
unchanged. The earlier overlapping partial-original cohort is excluded.

The output will reference existing qualified arrays under
`results/structural_alphabet/audited-esmfold-all-completed-20260922-v1/`.
No arrays are averaged or re-predicted. Resources are one CPU equivalent,
32 GiB RAM, no swap, 2 GiB output and a 50 GiB disk gate, with an uncalibrated
0.5–12 hour allowance after predecessors complete. Being queued is not a
completed qualification or evolutionary result. Paired-input reconstruction
and updated supported tree fits remain subsequent stages.
