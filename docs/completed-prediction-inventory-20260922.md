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


## Full refreshed AlphaFold residue readback queued

`scripts/readback_refreshed_afdb_residues.py` waits for the exact live mapping
controller and requires its successful catalog-agreement receipt before
reading the full mapped output. The service is
`fungal-current-afdb-residue-readback-20260922.service`; the pinned plan is
`metadata/current_afdb_residue_readback_plan.json`. It inherits the preparation
plan's complete frozen-input checks and writes to
`results/structural_markers/gdm-current-residue-readback-20260922-v1/`.

For every mapped marker/taxon record, the auditor checks the source amino-acid
sequence against its hash, the ungapped Stockholm alignment and the complete
mmCIF polymer. It requires exactly one Cα coordinate per sequence position and
valid pLDDT values. A cumulative non-gap projection independently reconstructs
the sequence position for each retained matrix column. Every exported residue
identity, matrix position, sequence position and confidence value must match,
with no missing or extra rows. All mapped models and marker links must be
covered, and each link's confidence count, mean and threshold fraction must
agree. The audit shares BioPython Stockholm/mmCIF parsers with the producer;
it does not independently assess alignment quality or model accuracy.

An end-to-end two-model synthetic test with seven hand-specified residue links
passed. Rehashed semantic corruptions were rejected for shifted positions,
altered confidence, missing rows, extra rows and duplicate Cα coordinates.
Reproduce with `python scripts/check_refreshed_afdb_residue_readback.py`.
The fixture evidence is recorded in
`metadata/current_afdb_residue_readback_fixture_checks.json`; this establishes
implementation checks, not a passed production audit.

The queued audit is capped at one CPU, 16 GiB RAM and zero swap, with a 50-GiB
free-disk gate and 1-GiB output allowance. The 0.5–12 hour interval after mapping
completion is an uncalibrated planning allowance. PAE qualification, native
structural features and evolutionary integration remain separate stages.


## Refreshed AlphaFold native coordinate features queued

The full 30,588-model, 15,949,929-residue catalog is queued for native 3Di
extraction, independent export identity checks and coordinate-feature
reconstruction. `scripts/advance_refreshed_afdb_features.py` waits for the exact
residue auditor process and requires its passing full-output receipt, matching
plan hash, mapping receipt and model count. It uses the existing extraction and
audit scripts unchanged; its dependency gate is specific to the refreshed
AlphaFold residue audit.

The service is `fungal-current-afdb-native-features-20260922.service`, with
`metadata/current_afdb_native_feature_plan.json` and launch evidence in
`metadata/current_afdb_native_feature_launch.json`. The controller output is
`results/structural_alphabet/current-afdb-feature-controller-20260922-v1/`;
native exports and audited coordinate encodings will be in
`results/structural_alphabet/native-gdm-current-20260922-v1/` and
`results/structural_alphabet/coordinate-gdm-current-20260922-v1/`.

The Foldseek executable, commit and four source-file hashes match the previously
qualified build. The job is capped at four CPU equivalents, four native threads,
16 GiB RAM and zero swap, with a 100-GiB output allowance and a 150-GiB free-disk
gate. The 1–24 hour interval after residue readback is uncalibrated planning.
No GPU inference or paid resources are used. The final stage verifies native
partner selection and coordinate descriptors and records full feature-context
pLDDT. PAE binding/qualification and evolutionary integration are still required;
this coordinate-only stage cannot by itself certify PAE-qualified states.


## Mapping-bound refreshed AlphaFold confidence queued

The final confidence controller waits for both the exact native-coordinate
controller and the exact PAE-prefetch controller. It requires their successful
receipts, matching plans, complete model counts and source bindings. It then:

1. Revalidates all cached PAE matrices into a receipt bound to the completed
   residue mapping, using the existing `retrieve_marker_pae.py`.
2. Runs the unchanged `qualify_native_pae.py` on all audited coordinate arrays.
3. Independently reads back every qualified six-residue PAE context using
   `scripts/readback_afdb_pae_qualification.py`.

The readback verifies compressed and uncompressed JSON checksums and model,
version, sequence-hash, length and URL provenance. It checks all pre-existing
coordinate arrays are unchanged, then accumulates the maximum across all 36
ordered pairs in each valid context without using the qualifier helper. Invalid
states must retain NaN PAE values. Every per-model and aggregate confidence
count must agree. Because AFDB PAE JSON does not itself embed the amino-acid
sequence, sequence identity is bound through the retrieval receipt and model
provenance; this is not direct sequence validation from the PAE payload.

An end-to-end synthetic six-residue asymmetric matrix passes with the expected
one confidence-qualified state. Corruptions of context maxima, invalid-state
sentinels, coordinate arrays, sequence provenance and uncompressed checksums
are rejected. Run `python scripts/check_afdb_pae_readback.py`; recorded results
are in `metadata/current_afdb_pae_readback_fixture_checks.json`.

`fungal-current-afdb-confidence-20260922.service` runs
`scripts/advance_refreshed_afdb_confidence.py` with
`metadata/current_afdb_confidence_plan.json`; its verified waiting process is
recorded in `metadata/current_afdb_confidence_launch.json`. The controller is
under `results/structural_alphabet/current-afdb-confidence-controller-20260922-v1/`.
Mapping-bound PAE will be under
`results/structural_pae/gdm-current-mapping-bound-20260922-v1/`; qualified arrays
will be under `results/structural_alphabet/audited-gdm-current-20260922-v1/`.

Resources are capped at two CPU equivalents, 16 GiB RAM and zero swap, with a
150-GiB free-disk gate and 100-GiB output allowance. The uncalibrated execution
allowance is 1–48 hours after both predecessors complete, not an ETA. This
stage uses existing predictions and local cached confidence data; no GPU
prediction or paid infrastructure is launched. Production qualification and
full-context readback are pending; their queueing does not establish calibrated
confidence, biological accuracy or completed evolutionary analyses.
# Refreshed AlphaFold paired alignments queued

The extended ESMFold cohort has completed its mapping-bound PAE export:
4,363 requested and verified models, zero failed exports, 27,440,436,824
compressed bytes. The manifest hash and mapping binding were checked and
the receipt archived as `metadata/esmfold_extended_full_pae_export_receipt.json`.
Its controller has advanced to confidence qualification. Original-context
readback and full-cohort encoding integration are still pending; export
completion does not establish confidence-qualified coverage.

The subsequent fitting controller is also queued as
`fungal-current-afdb-paired-fits-20260922.service`; its plan is
`metadata/current_afdb_paired_fits_plan.json`. It requires the exact paired
input controller to finish successfully and checks the completed input and
independent array-readback receipts before using any new alignment.
Every eligible marker receives an AA LG+F+G4 tree search with 1,000 SH-aLRT
and 1,000 UFBoot replicates with BNNI, then AF+G4, AF+F+G4 and LLM+G4
structural branch fits on that marker's fixed AA topology. This is at most
125 markers and 500 fits; actual dimensions are recorded after qualification.

The fitting producer, report auditor, structural models and IQ-TREE executable
are unchanged and hash-pinned. The controller allows four concurrent
single-thread fits, 2 GiB per fit, with a 16 GiB total memory cap, no swap,
and four CPU equivalents. Planning reserves 50 GiB output and requires
100 GiB free disk and 32 GiB available memory before execution. The
uncalibrated 1–336 hour runtime range is inherited from the full-cohort
planning allowance and is not a measured ETA. No GPU or paid resources are
used. The final report audit checks the fit grid, provenance, branch tables
and numerical warnings; it does not independently recompute likelihoods.
These are conditional point estimates. Branch uncertainty, model adequacy,
direct structural comparisons and ecological inference remain separate work.

The next CPU stage is registered as
`fungal-current-afdb-paired-inputs-20260922.service`, with plan
`metadata/current_afdb_paired_inputs_plan.json` and controller
`scripts/advance_refreshed_afdb_paired_inputs.py`. It waits for the exact
confidence-controller process identity and requires successful full-cohort
confidence qualification and context readback for all 30,588 models, bound
to the refreshed mapping and pinned confidence plan. A completed download
alone does not unlock preparation.

The unchanged `prepare_paired_phylogenetic_inputs.py` builds paired AA/3Di
inputs at `results/phylogeny/paired-inputs-gdm-current-20260922-v1`.
`readback_paired_inputs_from_encodings.py` then independently reconstructs
the emitted characters from qualified arrays. Final acceptance checks the
source bindings, all paired characters, unchanged masks and eligibility
rules, the same sequence matrix, and the complete marker universe. The
controller writes per-marker coverage changes against the older AlphaFold
alignments. Existing outputs remain available for sensitivity comparisons.

The stage has one CPU equivalent, 32 GiB RAM, no swap, no GPU, 10 GiB output
allowance and a 100 GiB free-disk gate. The uncalibrated runtime allowance is
0.5–24 hours after dependencies complete; it is not a measured ETA. The plan
pins 388 existing inputs and implementation files. Preparation and readback
are queued, not completed; fitted trees and ecological effects are separate
downstream work. The new Hydnum models are absent from the older qualified
snapshot, so current model-availability counts cannot stand in for completed
paired alignments.
