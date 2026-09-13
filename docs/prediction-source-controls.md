# Prediction-source controls and expanded mapping

Prediction pipelines are now explicit strata in marker mapping. The production default is provider **GDM**, tool **AlphaFold Monomer v2.0 pipeline**. Only exact provider/tool matches are eligible; there is no fallback to another pipeline when that source is missing. The deterministic confidence/version/model-ID ranking is applied within that source, and does not establish accuracy. The previous immutable structural snapshots remain historical checkpoints.

`map_marker_structures.py` freezes complete inventory records before mapping and audits source availability for every marker-protein link. A concurrent unfinished final record is excluded; malformed complete records raise an error. The frozen inventory is saved outside Git with a checksum. Alternatives remain in the inventory and are counted separately from source-eligible models. All 125 markers and the 526-taxon design remain in scope.

The expanded GDM mapping was launched with:

```bash
python scripts/map_marker_structures.py --output results/structural_markers/gdm-expanded-v1
```

Its terminal receipt establishes the final coverage; an existing output directory or partially written residue table does not establish completion. The previous 422-model snapshot continues to anchor the already executed comparisons until a new mapping, PAE and native-state audit are complete. An explicit provider-selection test verifies that a higher-confidence model from another tool cannot replace the requested source.

## Available same-sequence comparison

Every cross-pipeline, exact-sequence model pair in the frozen inventory was compared. At this checkpoint, there is **one distinct sequence** with such a pair: UniProt B9W6V0, 173 residues, GDM model AF-B9W6V0-F1 v6 and contributed ATBC model AF-0000000007270089 v1. The latter reports **ColabFold/AlphaFold2 Monomer v1.5.5 Pipeline**. This is the available subset of the full inventory, not a separately selected pilot or a representative sample of pipeline effects.

Full polymer identity, coordinate checksums, unambiguous complete CA positions, finite coordinates and confidence in [0,100] were verified. The literal polymer-sequence field is accepted when the canonical field is absent, with exact sequence equality and no guessing of modification codes. Both model-specific PAE matrices were retrieved and validated. The contributed PAE response is gzip-encoded; the reader now decodes transport compression before parsing and records the response checksum/encoding separately from the locally compressed JSON artifact.

At pLDDT ≥70 in both models, 171/173 residues remain. Their proper-rotation CA RMSD is **1.086 Å**, and mean PAE-filtered local CA-distance change is **0.193 Å**, retaining 3,331/3,337 eligible local pairs. At pLDDT ≥90, 141 residues remain, RMSD is 0.979 Å and local mean change is 0.163 Å. Confidence filtering changes the residue set; these differences do not show which model is more accurate.

Both pipelines use sequence-derived AlphaFold-family predictions. This comparison does not provide experimental validation, eliminate circularity, quantify systematic source effects across fungi, or establish equivalence of confidence scales. Broader cross-pipeline coverage, ESMFold controls and experimental structures remain required. In particular, do not treat pipeline-induced differences as evolutionary changes between identical sequences.

```bash
python scripts/compare_prediction_sources.py \
  --inventory results/structural_markers/gdm-expanded-v1/source_inventory.jsonl \
  --output results/structure_sources/available-pairs-v1
python -m unittest discover -s tests -p test_structure_source_selection.py
```

The comparison uses all available cross-pipeline pairs, records pLDDT thresholds 0/70/90, and separately labels insufficient coverage. Local CA pairs are within 15 Å in either model and separated by at least three sequence positions; directional PAE must be ≤10 Å in both models for the filtered metric. These are the same geometric definitions used elsewhere in the project. The one-pair resource estimate was recorded before execution: one CPU worker, cached/serial PAE retrieval, less than 1 GB RAM and 10 MB output expected, seconds to minutes, no paid resources.

Receipts and results are in `metadata/prediction_source_*`. The large frozen inventory and model/PAE data remain outside Git. A figure is not needed for this single available sequence; the complete three-row table is retained for reproducibility.

## Expanded model catalog and full-design coverage

The frozen inventory for the expanded mapper now has a separately validated model catalog: **13,153 distinct GDM models**, linked by complete sequence identity to **13,654 marker proteins in 322 taxa**. Every selected coordinate checksum was verified. This makes confidence-matrix retrieval independent of the slower alignment-to-coordinate mapping. Retrieval is running with two workers under the pre-launch estimate in `metadata/expanded_pae_resource_plan.json`; the conservative uncompressed JSON allowance is 56.8 GB. The 412 candidate cache receipts require readback before reuse.

The catalog and residue mapping have different receipts. Before prefetched confidence is used with the final mapping, their model identities, sequence identities, versions and coordinate hashes must agree. Once both producers finish, rerun the PAE retriever against the final mapping in a new output directory: its validated shared cache can supply the matrices while its new receipt binds them to that mapping. Catalog completion does not establish completed residue mapping, confidence-qualified coverage or native 3Di validation.

```bash
python scripts/prepare_marker_model_catalog.py \
  --inventory results/structural_markers/gdm-expanded-v1/source_inventory.jsonl \
  --output results/structural_markers/gdm-prefetch-catalog-v1
python scripts/retrieve_marker_pae.py \
  --snapshot results/structural_markers/gdm-prefetch-catalog-v1 \
  --output results/structural_pae/gdm-prefetch-v1
python scripts/audit_marker_model_coverage.py \
  --catalog results/structural_markers/gdm-prefetch-catalog-v1 \
  --inventory-inputs data/prediction_inputs/markers-full-inventory-v2 \
  --output results/structural_markers/gdm-coverage-audit-v2
```

The coverage audit includes all **526 taxa × 125 markers = 65,750 slots**, partitioned into four mutually exclusive states:

| State | Taxon–marker slots |
| --- | ---: |
| Downloaded GDM model in the frozen catalog | 13,654 |
| Reuse candidate outside this catalog | 17,764 |
| No candidate in the completed query snapshot | 28,422 |
| Marker sequence unavailable | 5,910 |

There are linked models for **304/501 fungal entries** and **18/25 outgroups**. All 204 taxa without catalog models remain in the tables. Aphelidiomycota, Calcarisporiellomycota, Sanchytriomycota and the Corallochytrea outgroup bin have no catalog models. Their recovered marker sequences have no candidate in the completed query snapshot; this is not evidence that the proteins or folds are biologically absent.

![Full-sampling marker model coverage](figures/expanded_marker_model_coverage.svg)

The gold category includes nominated retrievals, later downloaded models and potentially other pipelines outside this particular source catalog. Consequently, the prominent Ascomycota coverage advantage cannot be read as a biological difference: acquisition order and separately frozen checkpoints contribute to the pattern. Query scope and database coverage also limit the no-candidate category. Local predictions are not counted as downloaded GDM models. Sequence-identical model reuse preserves taxon links but does not create independent structural observations.

Each group uses an equal 125-marker denominator per taxon. Manifest lineage bins are descriptive and need not have equal taxonomic rank. Counts here are taxon–marker links, whereas the prediction-input queue counts unique sequences; those denominators must not be interchanged. Model availability also precedes the confidence, correspondence and alignment masks used in evolutionary inference.

The audit checks the exact marker universe, per-link protein/sequence identity, completed-inventory states, mutually exclusive count partition and every upstream receipt artifact. All output and script hashes passed independent readback, and the plotted figure was visually inspected. Versioned tables and receipts are `metadata/expanded_marker_*`; full catalog links and coordinates remain outside Git.

## Completed expanded residue mapping

The expanded mapper completed successfully with **4,895,692 matrix-to-structure residue links**, **13,654 marker proteins**, **13,153 models**, and **322 linked taxa**. There are at least four model-linked taxa for **124 markers**, before joint site/confidence filters. Completion and per-marker counts are recorded in `metadata/expanded_marker_mapping_receipt.json` and `metadata/expanded_marker_mapping_coverage.tsv`.

`verify_marker_catalog_mapping.py` checked every artifact in both receipts and exact agreement of the model/version identities, coordinate paths and checksums, complete-sequence hashes, source policy and all taxon–marker links. Its passing receipt is `metadata/expanded_marker_catalog_mapping_agreement.json`. The model-provenance JSON files have different row orders, so semantic identity is checked rather than assuming byte equality.

```bash
python scripts/verify_marker_catalog_mapping.py \
  --catalog results/structural_markers/gdm-prefetch-catalog-v1 \
  --mapping results/structural_markers/gdm-expanded-v1 \
  --output results/structural_markers/gdm-catalog-mapping-agreement-v1.json
python scripts/extract_structural_alphabet.py \
  --snapshot results/structural_markers/gdm-expanded-v1 \
  --output results/structural_alphabet/native-gdm-expanded-v1
```

Native Foldseek extraction was launched for all 13,153 models (6,910,765 full-model residues), using the previously pinned binary/source version and four threads. The prelaunch plan reserves 16 GB RAM and 5 GB additional output space on the existing host; original coordinates are reused through symlinks. See `metadata/expanded_native_3di_resource_plan.json`. Native extraction completion and native feature/confidence validation remain distinct gates. The prefetched PAE still requires completion and a validated receipt bound to this final mapping before downstream use.

Native extraction subsequently finished successfully. Independent export readback checked all **13,153 models**, **6,910,765 residues** and **69,107,650 descriptor values**: exact complete-sequence hashes, state alphabet/length, agreement between database FASTA and descriptor exports, finite ten-value feature vectors, and complete model coverage. There are 26,306 zero-log-sequence-offset rows; their exported letters are not accepted as valid structural observations. The native output occupies approximately 866 MB, within the planned disk envelope.

```bash
python scripts/audit_native_exports.py \
  --native results/structural_alphabet/native-gdm-expanded-v1 \
  --snapshot results/structural_markers/gdm-expanded-v1 \
  --output results/structural_alphabet/native-gdm-expanded-exports-v1.json
```

`metadata/expanded_native_export_audit_receipt.json` pins the checked exports; `metadata/expanded_native_3di_config.json` records commands and software/source identities. This completed export audit precedes coordinate-based reconstruction of native features/partners and the six-residue pLDDT/PAE audit, which remains pending completed confidence acquisition. The existing small-snapshot evolutionary benchmarks have not been rerun on this expansion yet.

## Coordinate audit independent of confidence acquisition

`audit_3di_features.py` now supports a coordinate-only stage by omitting `--pae`. It reconstructs the spatial partner and all ten native features from the coordinates, validates the native validity mask and descriptor tolerance, and records focal and six-residue minimum pLDDT. Its receipt has status `complete_native_3di_coordinate_audit`; it has no PAE receipt or joint-confidence count, and the NPZ files omit `feature_max_pae`. Missing confidence is not represented as zero or accepted as passing.

`qualify_native_pae.py` then joins a completed coordinate audit to completed PAE acquisition only when their mapping receipt hashes and complete model/version universes agree. It checks every encoding and PAE artifact hash and complete-sequence identity, then computes the maximum over all 36 directional entries in each six-residue feature context. Invalid sites remain NaN. It produces the same final encoding fields and confidence counts as the original integrated audit.

The full historical 422-model dataset was reprocessed through both stages. `compare_native_audits.py` verified exact agreement of **all 2,954 arrays**, every model summary and all totals with the earlier completed audit, including NaN placement. The reference contains 215,518 residues and 149,351 valid states passing both six-residue pLDDT70 and PAE10. Three geometry tests and two additional directional-context tests passed. This is a computational regression across an existing full checkpoint, not a new biological pilot or evidence of expanded-audit completion.

```bash
python scripts/audit_3di_features.py \
  --native results/structural_alphabet/native-marker-v2 \
  --snapshot results/structural_markers/snapshot-v2 \
  --output results/structural_alphabet/coordinate-regression-v1
python scripts/qualify_native_pae.py \
  --coordinates results/structural_alphabet/coordinate-regression-v1 \
  --pae results/structural_pae/snapshot-v1 \
  --output results/structural_alphabet/qualified-regression-v1
python scripts/compare_native_audits.py \
  --reference results/structural_alphabet/audited-marker-v2 \
  --candidate results/structural_alphabet/qualified-regression-v1 \
  --output results/structural_alphabet/staged-native-regression-v1.json
```

After this successful regression and a recorded resource estimate, the expanded coordinate-only audit was launched with one worker and numerical-library threads fixed at one:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python scripts/audit_3di_features.py \
  --native results/structural_alphabet/native-gdm-expanded-v1 \
  --snapshot results/structural_markers/gdm-expanded-v1 \
  --output results/structural_alphabet/coordinate-gdm-expanded-v1
```

Its 13,153 models contain 6,910,765 residues; maximum model length is 2,413. The plan allows 8 GB RAM, 2 GB additional output and approximately 0.5–3 hours on the existing host, with no new charges. Coordinate parsing and quadratic within-model partner calculations can now overlap PAE retrieval. The existing output directory is immutable; a partial directory is not a completed audit or an automatically resumable result. Final mapping-bound PAE acquisition and both completed audit stages are still required before expanded phylogenetic inference. Logs: `logs/native_coordinate_gdm_expanded_v1.log`.

## Expanded coordinate reconstruction completed

The coordinate-only production audit finished successfully for **all 13,153 models** and **6,910,765 residues**. It reproduced the pinned native partner assignments and all ten features within the declared four-significant-digit descriptor tolerance. There are **6,884,459 valid states** and **26,306 invalid terminal states**. Of the valid states, 5,136,862 meet focal pLDDT70 and **4,719,638** meet the minimum pLDDT70 across the complete six-residue feature context.

Independent readback verified every one of the 13,153 encoding checksums, unique model/version counts, all summary totals and the producer script hash. The output occupies approximately 95 MB. The completion receipt is `metadata/expanded_coordinate_audit_receipt.json`; encodings and the per-model summary remain under `results/structural_alphabet/coordinate-gdm-expanded-v1`.

These counts cover full-model residues, not just retained alignment sites. The receipt explicitly has no PAE parent, the coordinate-only NPZs omit `feature_max_pae`, and no joint pLDDT/PAE count is reported. Final mapping-bound PAE acquisition and `qualify_native_pae.py` remain required before rebuilding expanded matched AA/3Di phylogenetic inputs. Passing coordinate reconstruction is not experimental validation or an evolutionary acceleration result.

## Dedicated AlphaFold–ESMFold controls prepared

The first 4,895 completed local ESMFold receipts had no exact-sequence overlap
with the frozen expanded GDM mapping: the production queue targets missing
structures. A dedicated control input now contains 266 unique full proteins
(78,593 residues), 270 marker/taxon links from 126 taxa, and 23 represented
manifest lineage/role groups. All associated source identities are retained.

Selection takes up to four sequences per study-role/manifest-lineage/128-residue
length-bin/AlphaFold-mean-pLDDT-bin stratum. The confidence bins are below 70 and
at least 70. Candidates must be canonical full sequences of at most 512 residues;
salted SHA256 ranking makes selection deterministic. The union deduplicates
sequences, while the selection table retains their potentially multiple strata.
This is an availability- and length-limited predictor control, not a random
sample of fungal proteins or an independent replicate for every taxon link.

Reproduce with `python scripts/prepare_predictor_control_inputs.py --output
 data/prediction_inputs/predictor-controls-v1` using a new output location if
already present. Reference model metadata, candidate FASTA and full links remain
in that input directory; compact receipts and selection/link tables are
versioned under `metadata/predictor_control_*`. Readback checked every output
hash, canonical sequence/hash identity, exact reference sequence universe and
all 266 reference coordinate checksums.

The control predictions have not launched. Based on 4,919 completed production
receipts, length-bin median inference times imply 0.467 GPU hours; a twofold
allowance plus ten minutes of startup gives a 1.10-hour planning envelope,
with 24 GB GPU memory and 2 GB disk headroom. This is an extrapolation, not a
completion guarantee. After the active missing-model production finishes,
recheck GPU occupancy and use the same pinned predictor in a separate output
directory. Exact-sequence geometry and confidence comparisons remain pending.
Both methods derive structures from sequence; agreement alone cannot resolve
shared training biases or replace experimental validation.

### Control predictions launched

A fresh device check found GPU 0 idle with 41 MiB used and no compute processes;
the earlier unrelated workload had ended. Launched all 266 controls on GPU UUID
`GPU-56e78ad3-b4d7-54f6-38fa-e95729009960`, while missing-model production
continues on GPU 1. Configuration differs from the existing prediction run
only in the input receipt and GPU UUID. Both cards have the same hardware model.
The first four control predictions completed successfully; full completion and
independent artifact readback remain pending. Current launch resources and
configuration are pinned in `metadata/predictor_control_launch_resource_plan.json`
and `metadata/predictor_control_prediction_config.json`.

```bash
CUDA_VISIBLE_DEVICES=GPU-56e78ad3-b4d7-54f6-38fa-e95729009960 \
/home/bizon/anaconda3/envs/esmfold2/bin/python scripts/run_marker_predictions.py \
  --inputs data/prediction_inputs/predictor-controls-v1 \
  --checkpoint data/prediction_models/esmfold-v1 \
  --output results/predictions/esmfold-controls-v1 --max-length 512 --limit 266
```

After the full control chunk completes, run the independent auditor with
`--predictions results/predictions/esmfold-controls-v1`,
`--inputs data/prediction_inputs/predictor-controls-v1`,
`--links data/prediction_inputs/predictor-controls-v1/reference_links.tsv`, and a
new `--output` directory. The auditor accepts either full sequence IDs or
sequence hashes, rejects inconsistent identities, and records the link-table
hash. Three focused tests cover both formats and malformed/conflicting links.
This adds input compatibility to the existing artifact audit; the tests do not
establish accuracy of future predictions.

### Exact-sequence geometry comparison prepared

`compare_esmfold_controls.py` requires the complete control chunk and its
independent artifact audit before analyzing every reference pair. It verifies
full residue identities/numbering, coordinate and prediction hashes, and
mapping-bound reference PAE. It evaluates jointly retained pLDDT thresholds
0, 70 and 90; fewer than 50 residues or less than half of the full protein
produces an explicit coverage exclusion with missing metrics.

For accepted comparisons, proper-rotation CA superposition gives RMSD. The
local metric compares residue-pair distances within 15 Å in either prediction,
excluding pairs separated by fewer than three sequence positions. Additional
summaries require PAE at most 5, 10 or 15 Å in both directions in both models.
These filters change the residue/pair sets; differences across thresholds do
not by themselves diagnose prediction error. This local metric is not lDDT,
and predictor differences are not biological branch lengths.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/compare_esmfold_controls.py \
  --inputs data/prediction_inputs/predictor-controls-v1 \
  --predictions results/predictions/esmfold-controls-v1 \
  --audit results/predictions/audit-controls-v1 \
  --pae results/structural_pae/gdm-expanded-v1 \
  --output results/prediction_controls/esmfold-af-comparisons-v1
```

This command is prepared, not yet executed. The planned comparison uses one CPU
worker, up to 4 GB memory and 0.1 GB output headroom. Four geometry/filter tests
pass: rigid transformations produce zero change, reverse-direction PAE in
either model excludes the affected pair, insufficient coverage stays missing,
and nonfinite confidence is rejected. Full empirical comparison, sensitivity
interpretation and experimental validation remain pending.

## Full 266-control prediction comparison completed

All 266 control predictions finished without interruption, OOM or deferred
entries. Independent readback passed every model, 270 source links, 126 taxa
and 81 marker identities. The geometry comparison completed all 798
protein/threshold rows. The pLDDT ≥70 rule retains sufficient matched residues
for 209 controls and excludes 57; pLDDT ≥90 retains 88 and excludes 178.

| Joint pLDDT cutoff | Compared proteins | Median CA RMSD (Å) | Median local distance change (Å) | Median PAE≤10 local change (Å) |
| --- | ---: | ---: | ---: | ---: |
| 0 | 266 | 11.848 | 0.891 | 0.316 |
| 70 | 209 | 1.114 | 0.236 | 0.219 |
| 90 | 88 | 0.535 | 0.157 | 0.157 |

These rows change both the protein cohort and residue sets. Holding the 209
proteins that qualify at 70 fixed, their full-protein median RMSD is 7.615 Å,
compared with 1.114 Å after residue filtering. The median within-protein change
is −3.545 Å. Holding the 88 proteins qualifying at 90 fixed, the corresponding
full-protein and filtered medians are 3.840 and 0.535 Å. These comparisons still
change residue positions. They describe sensitivity of predictor agreement to
confidence selection, not improvement in accuracy, experimental validation,
biological change or independent-replicate statistical significance. Controls
were intentionally stratified by availability, lineage, length and AF confidence;
these medians are not estimates for all fungal proteins.

`compare_esmfold_controls.py` and `summarize_predictor_controls.py` produced the
completed artifacts under `results/prediction_controls/esmfold-af-comparisons-v1`
and `esmfold-af-summary-v1`. The summary checks the complete 266-by-three grid,
finite/nonnegative geometry and nested directional PAE retained-pair counts.
Receipts, all comparison rows and both changing-/same-cohort summaries are
versioned in `metadata/predictor_control_*`. Leading disagreements still require
domain/orientation review and checks against experimental structures.
