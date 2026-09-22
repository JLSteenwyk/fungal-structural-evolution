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

All five cohorts and their full combined inventory have completed residue
mapping. Independent combined readback verified 8,335,615 matrix-residue links,
25,509 marker/protein links and all 25,322 models. Native feature/confidence
qualification, regenerated paired alignments, updated phylogenetic fits and
evolutionary tests remain unfinished.

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
The mapping gate has passed. Confidence-qualified encoding integration remains
a separate gate before rebuilding paired phylogenetic inputs. Completion
receipts are archived under `metadata/esmfold_all_completed_mapping_completed.json`,
`metadata/esmfold_all_completed_residue_readback.json` and
`metadata/esmfold_all_completed_mapping_controller_receipt.json`.

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

## Availability across the completed ESMFold and frozen AlphaFold catalogs

The full marker universe contains 59,840 recovered marker/protein records
(58,883 distinct sequences), across 125 markers and 526 sampled entries.
An exact marker/taxon/protein/sequence join partitions every recovered record:

| Availability in these catalogs | Marker/protein records |
|---|---:|
| ESMFold only | 25,277 |
| AlphaFold only | 13,422 |
| Both | 232 |
| Neither | 20,909 |

Together, these catalogs cover 38,931 records (65.1%), with at least one
model-linked marker in every one of the 526 sampled entries. These are
pre-qualification counts; they do not supersede confidence-qualified coverage
counts or demonstrate adequate within-family or within-lineage sampling.
They also do not represent all proteins in the genomes. The denominator
excludes unrecovered marker sequences, which remain separately tabulated.

“Neither” means absent from the completed ESMFold inventory and the frozen
GDM/AlphaFold catalog used in this comparison. It does not establish absence
from newer cached downloads or public databases, and must not be converted
directly into a prediction queue without a current retrieval check. Models
shared by exact sequences are not independent observations.

Reproduce with `python scripts/summarize_completed_marker_availability.py
--plan metadata/completed_marker_availability_plan.json` after selecting a fresh
output directory in a copied plan. The plan pins all inputs and the producer,
and records a one-CPU, 1-GiB memory, 0.1-GiB output, 1–5-minute planning
allowance. This computation uses no GPU. Every source link is checked against
the global inventory and its source receipt. A separate pandas two-left-join
calculation reproduced all 59,840 output classifications, all 526 taxon rows and
the aggregate counts; it also checked marker/taxon uniqueness.

The receipt and per-taxon table are versioned as
`metadata/completed_marker_availability_receipt.json` and
`metadata/completed_marker_taxon_availability.tsv`. The complete per-marker
availability table remains in
`results/structures/completed-marker-availability-20260922-v1/`.


## Original-cohort confidence processing completed

The complete 10,522-model original cohort has passed native feature extraction,
coordinate reconstruction, confidence qualification and independent checks
against its original directional PAE arrays. Across 3,355,507 residues,
3,334,463 have valid structural-alphabet states. Of these, 2,602,366 pass the
focal-residue pLDDT threshold of 70; 2,427,934 pass that threshold over the full
six-residue feature context; and 2,424,268 also pass the directional-context PAE
threshold of 10 Å. These counts are per model residue, before marker-alignment
filtering, and are not independent observations or measures of prediction
accuracy.

The PAE readback checked all 36 ordered pairs within each valid six-residue
context. All controller stage hashes and every model-summary count were
additionally checked for consistent totals and nested masks. Archived receipts:

- `metadata/esmfold_original_full_feature_completed_receipt.json`
- `metadata/esmfold_original_full_qualified_receipt.json`
- `metadata/esmfold_original_full_pae_context_readback.json`

Long-cohort confidence processing has also completed; extended-cohort PAE export remains active.
The existing qualified-union controller will advance only after their full
qualification and readback gates pass.

## Refreshing downloaded AlphaFold marker availability

The retrieval process has continued beyond the September 13 frozen catalog.
A new immutable snapshot captures 1,316,471 complete retrieval-log records
(1,542,870,552 bytes), including all statuses and repeated accessions. These
records are not counts of unique accessions or available models. The snapshot
is `results/structures/afdb-inventory-20260922-v1/inventory.jsonl`; its hash and
source-prefix readback are archived in
`metadata/current_afdb_inventory_freeze_receipt.json`.

`scripts/freeze_append_only_inventory.py` copies only the complete-line prefix
within the source's initial byte length, then independently rereads that source
prefix. Concurrent later appends are excluded. The complete-line and unfinished-
tail cases passed focused checks. The unchanged
`scripts/prepare_marker_model_catalog.py` completed on this frozen input under
`fungal-current-afdb-marker-catalog-20260922.service`, selecting models with the
same GDM/AlphaFold Monomer v2.0 source policy as the earlier catalog and verifying
selected coordinate hashes. Output is
`results/structural_markers/gdm-current-catalog-20260922-v1/`.

The resource plan is `metadata/current_afdb_marker_catalog_plan.json`: one CPU,
32 GiB RAM, no swap, 3 GiB output allowance and a 50-GiB free-disk gate. The
0.1–2 hour interval is uncalibrated planning. The launch record is
`metadata/current_afdb_marker_catalog_launch.json`. No new predictions or paid
resources are involved. Updated combined coverage and the actual remaining
prediction queue require the completed catalog and independent joins; the
older 20,909-record gap remains a statement about the older frozen catalogs.


## Current combined availability: 93.9% of recovered marker records

The refreshed GDM catalog contains 30,588 selected models linked to 31,384
marker/protein records across 339 taxa. Its coordinate hashes passed catalog
construction checks. `scripts/readback_marker_catalog_selection.py` independently
replayed all 1,316,471 frozen log records and reproduced every selected model and
marker link, including latest-accession status, source policy and tie handling.
The readback does not repeat coordinate parsing or assess biological accuracy.

Together with the 25,322 completed ESMFold models, the updated catalogs cover
**56,171 of 59,840 recovered marker/protein records (93.87%)**:

| Current catalog availability | Marker/protein records |
| --- | ---: |
| ESMFold only | 24,787 |
| AlphaFold only | 30,662 |
| Both | 722 |
| Neither | 3,669 |

All 526 sampled entries have at least one model-linked marker. The refreshed
AlphaFold catalog retains every earlier AlphaFold-linked marker record and
adds 17,730 links. Of those additions, 17,240 fill previously uncovered records;
the other 490 overlap ESMFold coverage. The earlier 65.1% figure remains valid
for its older frozen catalog but is superseded for current downloaded-model
availability. Neither percentage is proteome-wide or confidence-qualified
coverage, and existing evolutionary analyses have not yet incorporated these
additional models. The 3,669 uncovered records still require checks against
other source policies and later retrievals before defining a prediction queue.

The unchanged `summarize_completed_marker_availability.py` ran under
`metadata/current_marker_availability_plan.json`. An independent pandas join
reproduced all 59,840 classifications and all 526 taxon summaries, and confirmed
that no earlier AlphaFold marker link was lost. Current outputs are under
`results/structures/current-marker-availability-20260922-v2/`. Archived evidence:

- `metadata/current_afdb_marker_catalog_receipt.json`
- `metadata/current_afdb_marker_catalog_selection_readback.json`
- `metadata/current_marker_availability_receipt.json`
- `metadata/current_marker_availability_join_readback.json`
- `metadata/current_marker_taxon_availability.tsv`

The updated AlphaFold models require residue mapping, confidence processing
and integration before downstream phylogenetic and structural comparisons.


## Refreshed AlphaFold mapping and PAE preparation running

Two independent stages now run through
`scripts/advance_refreshed_afdb_inputs.py`, with 392 input/script pins recorded
in `metadata/current_afdb_input_preparation_plan.json`:

- `fungal-current-afdb-mapping-20260922.service` maps all 30,588 selected models
  and 31,384 marker/protein links onto the complete frozen profile matrix. It
  then requires exact catalog/mapping identity agreement. Output is
  `results/structural_markers/gdm-current-mapping-20260922-v1/`.
- `fungal-current-afdb-pae-20260922.service` retrieves or revalidates all
  version-matched PAE matrices. Output is
  `results/structural_pae/gdm-current-prefetch-20260922-v1/`. A failed model
  remains explicit in the manifest and prevents acceptance of the full stage.

The catalog contains 15,949,929 protein residues, a maximum length of 2,416,
and 10,859,772,095 PAE matrix entries. At planning time 13,127 cached PAE
receipt files were present and 17,461 were absent. File presence is not
validation: every reused matrix is checked for provenance, compressed and
uncompressed checksums, dimensions and valid numeric entries.

Mapping is capped at one CPU and 64 GiB RAM; PAE retrieval at two CPU
equivalents and 8 GiB RAM, with two retrieval threads. Both have zero swap.
The plan reserves 30 GiB for mapping and 200 GiB for PAE cache/output, gated on
300 GiB free disk. Uncalibrated runtime allowances are 0.5–12 hours for mapping
and 2–48 hours for PAE preparation; these are not measured ETAs. Existing local
resources are used without new GPU inference or paid infrastructure. Launch
identities are recorded in `metadata/current_afdb_input_preparation_launch.json`.

Independent residue readback, mapping-bound PAE checks, native structural
features and confidence qualification remain subsequent requirements. Passing
these preparation stages alone will not establish new evolutionary results.

## Long ESMFold cohort confidence completed

All 5,510 models of length 513–768 residues now pass the feature/confidence
pipeline and original directional-PAE-context readback. Of 3,386,360 model
residues, 3,375,340 have valid structural-alphabet states, 2,069,197 pass pLDDT
70 over the full feature context, and 2,065,046 additionally pass maximum
context PAE 10 Å. All stage hashes and per-model nested-mask counts were checked
before archival. Counts precede marker-alignment filtering and do not measure
biological accuracy. Receipts are archived as
`metadata/esmfold_long_feature_completed_receipt.json`,
`metadata/esmfold_long_qualified_receipt.json` and
`metadata/esmfold_long_pae_context_readback.json`.
