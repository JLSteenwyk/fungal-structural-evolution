# Local structure prediction

The first local method is the Hugging Face implementation of ESMFold v1,
`facebook/esmfold_v1` revision `ba837a39b67e59941c3f017d6c2a064f567038d9`.
The official [model card](https://huggingface.co/facebook/esmfold_v1) documents
single-sequence prediction without an MSA or external database lookup.
`prepare_esmfold_checkpoint.py` verifies the cached safetensors SHA256 against
the publisher revision API and the small configuration files as Git blobs.
The tracked receipt contains pinned URLs to reconstruct the checkpoint.

`prepare_marker_prediction_inputs.py` freezes complete records from the ongoing
UniProt/AFDB inventories, validates marker inputs, and preserves all 59,840
marker-to-taxon links. Of 58,883 unique sequences in the first snapshot, 20,480
lack a candidate in completed taxon-specific cross-reference inventories,
13,170 have nominated retrievals pending, 8,464 have verified reuse receipts,
and 16,769 await inventory. These states do not prove absence from all
structural databases. Reuse receipts route the queue; downstream structural
use still requires validation of the actual coordinates.

The initial production chunk requests 128 complete proteins of at most 512
residues from a round-robin taxon queue, shorter proteins first within each
taxon. This is part of the full analysis, not a separate biological pilot.
The chunk is biased toward short proteins and cannot establish performance
for longer proteins. All longer/noncanonical proteins remain explicitly
deferred. New snapshots are required as inventories finish. The marker queue
is not the complete proteome atlas.

Before launch, the resource envelope is one existing 48 GB GPU, four CPU
threads, approximately 40 GB host RAM, and 2 GB output headroom for 128 proteins
up to 512 residues. The 8.44 GB checkpoint is already cached. A provisional
planning range of 5–300 seconds per short sequence means 11 minutes–11 hours
for this chunk, excluding initialization. This is not measured performance.
Peak allocation and inference time are recorded per protein. Long proteins
need a separate measured resource envelope. No paid infrastructure or remote
inference endpoint is used. Check GPU availability before each launch;
GPU 0 has an unrelated active workload.

The environment directory name is `esmfold2`, but the chosen checkpoint is
ESMFold v1. The installed package inventory is recorded in
`environments/esmfold-local-packages.json`; this inventory alone is not a
reconstructable environment lock. Exact model configuration, versions and
modeling source hash are recorded with each run.

```bash
CUDA_VISIBLE_DEVICES=GPU-21f34861-216e-5737-10cd-aad35dc83ac1 \
  /home/bizon/anaconda3/envs/esmfold2/bin/python scripts/run_marker_predictions.py \
  --inputs data/prediction_inputs/markers-v1 \
  --checkpoint data/prediction_models/esmfold-v1 \
  --output results/predictions/esmfold-marker-v1 --max-length 512 --limit 128
```

Settings: FP16 ESM stem, FP32 folding trunk/heads, attention chunks of 64,
one sequence per batch, checkpoint default four trunk passes, eval/inference
mode, seed 20260913 and disabled TF32. Bitwise reproducibility across platforms
is not claimed. Predictions are unrelaxed isolated chains without templates
or an MSA. Transformers native confidence is 0–1; the exporter multiplies by
100 for PDB B factors and CA confidence. Validation checks the complete
sequence, consecutive residue numbering, finite CA coordinates, confidence
range/rounding, and square finite nonnegative PAE. Full PAE and unrounded
CA confidence are saved in compressed NumPy files.

Per-protein checksums and configuration hashes permit validated resume.
An out-of-memory failure records a deferred protein and ends the chunk
without trimming its sequence. SIGTERM stops after the current prediction.
A chunk completion receipt does not establish completion of the full queue.
ESMFold predictions remain separate from AFDB acquisition; confidence measures
are not interchangeable calibrated errors. Source effects, experimental
controls and prediction circularity remain to be assessed.

## First production chunk and expansion

All 128 requested predictions finished, with 129 taxon-marker links across 129 taxa and 15 markers. Independent Biopython/PDB and NumPy readback passed. Protein lengths were 70–255 residues (median 119), inference times 0.845–2.899 seconds (median 0.944), and maximum measured GPU allocation 8.92 GB. Mean CA pLDDT ranged 38.81–91.81; low-confidence models remain labeled, not discarded or interpreted as validated folds. An execution audit checked that the only missing checkpoint parameters were the contact-regression bias/weight, installed a rejecting hook on that head, and completed folding inference without invoking it. Repeated PAE and CA confidence were identical for the checked sequence.

Expansion retains the same configuration and resumes all 10,394 remaining canonical candidates no longer than 512 residues. A provisional 1–30 seconds per sequence gives approximately 3–87 GPU hours plus serialization and initialization; the first short-protein measurements do not validate the upper length regime. Reserve 50 GB output headroom and the same single GPU/four CPU threads. Per-sequence memory/runtime receipts and an explicit OOM stop remain in force. This expansion still leaves 9,958 length/alphabet-deferred sequences and future inventory/whole-proteome candidates pending.

## Prepared follow-on missing-model inputs

A runner-compatible follow-on input snapshot now covers all **7,799 additional unique sequences** identified by the completed full-inventory delta. It preserves **7,895 taxon–marker links across 156 taxa**. Every sequence is disjoint from the original 20,480-candidate queue, including that queue's deferred long proteins. A newly frozen AFDB inventory and the current local prediction directory were checked for verified exact-sequence reuse; none of these additional candidates had a reusable model at this checkpoint. This is a query-snapshot statement, not global database absence.

The prepared dispositions are **4,252 short/canonical prediction candidates**, **3,495 canonical length-deferred candidates**, and **52 noncanonical candidates**. Noncanonical status takes precedence over length in this preparation table; total deferred proteins are 3,547. Complete sequences and links remain intact. The running GPU job continues with its original frozen inputs; no additional GPU run has been launched or queued by these preparation scripts.

```bash
python scripts/prepare_additional_prediction_inputs.py \
  --previous data/prediction_inputs/markers-v1 \
  --current data/prediction_inputs/markers-full-inventory-v2 \
  --delta data/prediction_inputs/full-inventory-delta-v2 \
  --existing-predictions results/predictions/esmfold-marker-v1 \
  --output data/prediction_inputs/markers-followon-v1
python scripts/estimate_followon_prediction_runtime.py \
  --inputs data/prediction_inputs/markers-followon-v1 \
  --predictions results/predictions/esmfold-marker-v1 \
  --output results/predictions/followon-resource-v1
```

The length-stratified resource projection freezes **3,640 observed prediction receipts**, using separate 1–128, 129–256, 257–384 and 385–512 residue bins. Applying each bin's median runtime to the actual follow-on lengths gives **8.23 GPU-hours** of inference. The analogous sums using observed 10th/90th-percentile times are 5.98/11.18 hours; these are planning scenarios, not confidence intervals. A 50% allowance above the latter gives **16.77 hours** for scheduling, including overhead and contention. The plan allows 24 GB VRAM and 20 GB output space on the existing GPU, without new charges. These estimates do not apply to deferred long proteins or the whole multi-million-protein atlas.

Before execution, confirm the current GPU run has finished, revalidate new reuse evidence and GPU availability, and use a new prediction output directory because the input receipt changes. No further user approval is required for continuing on the existing authorized host. The current runner's exact configuration and cache hashes remain protected. After prediction, independently validate artifacts and retain source strata before structural comparisons.

All prepared input hashes, sequence identifiers, disjointness and short-candidate counts passed independent readback. Receipts and runtime-bin summaries are `metadata/followon_prediction_*`; the complete source/taxon links and timing observations remain under their recorded data/result directories. No prediction-completion claim follows from these preparation receipts.

## Independent readback of an active production snapshot

The auditor now supports `--snapshot-live` to freeze the current list of
per-model JSON receipts and independently verify their source sequence,
PDB coordinates/numbering, NPZ confidence/PAE and artifact hashes while the
producer continues. Per-model receipts are published only after their artifacts
are complete. The resulting audit is explicitly a partial prediction snapshot;
it does not use an earlier `last_chunk.json` as evidence of current production
completion, and reports remaining eligible count and chunk-receipt hash as
unknown. The receipt TSV records every included prediction-receipt checksum.

Started a 5,121-receipt snapshot audit under
`results/predictions/audit-marker-partial-v2`, with output in
`logs/esmfold_partial_snapshot_v2.log`. It has not yet completed. All original
marker/taxon links are retained for included sequence identities. Following
readback, this snapshot can support a separate ESMFold residue-mapping and
structural-confidence pipeline; candidate availability alone does not establish
usable paired coverage or justify mixing predictor effects into branch rates.

```bash
python scripts/audit_local_predictions.py \
  --predictions results/predictions/esmfold-marker-v1 \
  --inputs data/prediction_inputs/markers-v1 \
  --output results/predictions/audit-marker-partial-v2 --snapshot-live
```

The existing completed-chunk audit remains the default. Two integration tests
use a stale historical chunk receipt to verify that partial snapshots remain
labeled partial and completed-chunk mode still rejects inconsistent counts.
Three link-identity tests also pass. One CPU worker, 4 GB memory and 0.1 GB
output headroom are planned; no new inference or paid resources are needed.

## Partial ESMFold readback completed; mmCIF conversion running

All 5,121 models in the frozen snapshot passed independent artifact validation:
5,161 marker/taxon links, 201 taxa and 86 markers. The included taxa comprise
122 Ascomycota, 76 Basidiomycota, two Aphelidiomycota and one
Basidiobolomycota. This is the composition of the current snapshot, not the full
queue or confidence-qualified structural coverage. In particular, Aphelidiomycota
model availability must still pass residue mapping and feature confidence
before claiming that its paired-coverage gap has been filled.

`convert_esmfold_snapshot.py` is now exporting this audited snapshot into
`results/structural_inventory/esmfold-partial-v1`. It writes explicit canonical
polymer sequences and preserves every atom field, residue ID, coordinate and
PDB confidence decimal in mmCIF, then parses each output back and compares all
fields. It rejects alternate locations, insertion codes, unexpected heteroatoms
and inconsistent sequence/numbering. Three focused conversion tests pass, and
over 500 actual models have converted successfully. Full conversion is pending.

```bash
python scripts/convert_esmfold_snapshot.py \
  --audit results/predictions/audit-marker-partial-v2 \
  --predictions results/predictions/esmfold-marker-v1 \
  --output results/structural_inventory/esmfold-partial-v1
```

The local inventory uses sequence-based record IDs, explicit ESMFold provenance,
configuration-derived model IDs and a local representation version. It retains
original PDB/NPZ paths and checksums. No UniProt accession or remote PAE URL is
invented. The mapper now accepts `--inventory` and local `record_id` fields;
its existing exact provider/tool restriction remains enforced. Three source/
identity tests pass. After conversion completes, a separate local ESMFold mapping
can use `--provider local --tool 'ESMFold v1'`; local PAE binding and native
feature qualification still require downstream work. Existing AlphaFold
snapshots remain unchanged.

## Local PAE export for structural-feature qualification

`export_local_marker_pae.py` is prepared for the completed local ESMFold residue
mapping. It requires a local/ESMFold-only source policy, validates every mapped
model's original prediction receipt and NPZ hash/configuration/sequence identity,
and exports directional PAE to gzip JSON in the downstream matrix schema.
Every exported array must agree exactly with the original NPZ after decompression
and parsing. No symmetrization, rounding, remote download or AlphaFold provenance
is introduced. Matrix dimensions, finite/nonnegative values and the original
ESMFold declared maximum are checked before export.

The output manifest explicitly labels local prediction exports and records the
original NPZ and prediction/configuration hashes. Its receipt binds the entire
model set to the exact residue-mapping receipt, allowing native feature
qualification to use the existing identity checks. Three tests cover directional
precision-preserving roundtrips, wrong sequence identity and invalid matrices or
maximum values. This stage has not launched: mmCIF conversion and completed
local mapping are prerequisites. The maximum current input envelope is 5,121
models and 303,017,975 matrix entries, with one CPU worker and conservative
resource estimates in `metadata/local_pae_export_resource_plan.json`.

## Follow-on missing-model predictions launched

After the 266 controls finished, GPU 0 was again independently verified free.
The prepared follow-on input has 4,252 canonical sequences no longer than 512
residues; none has a model in the freshly checked AlphaFold cache, and the input
is disjoint from the entire original prediction queue. Launched all 4,252 in
`results/predictions/esmfold-followon-v1` using the same pinned configuration;
only the input receipt differs from the completed control run. The first model
completed successfully. GPU 1 continues the original queue.

The existing empirical forecast is 8.23 GPU hours of inference, with a 16.77-hour
planning envelope and 20 GB output headroom. These are planning scenarios, not
completion guarantees. Fresh cache/device observations, source hashes and run
configuration are recorded in `metadata/followon_prediction_launch_*`. The 3,547
long/noncanonical input entries remain explicitly deferred by this configuration.
No paid resources were provisioned.

## All 5,121 mmCIF conversions completed; local mapping started

Conversion finished for every audited model. Independent readback verified the
completed inventory/provenance hashes, unique model IDs and every converted CIF
checksum. The producer's per-model roundtrip preserved all exported atom fields
and canonical sequences. The completion receipt is versioned as
`metadata/esmfold_conversion_receipt.json`.

Started the exact-sequence local ESMFold mapping across the full 526-taxon,
125-marker source design:

```bash
python scripts/map_marker_structures.py \
  --inventory results/structural_inventory/esmfold-partial-v1/inventory.jsonl \
  --provider local --tool 'ESMFold v1' \
  --output results/structural_markers/esmfold-partial-v1
```

The mapping process is live; its log is `logs/esmfold_marker_mapping_v1.log`.
This uses the same verified profile alignments and retained matrix columns as
the separate GDM mapping, with provider-specific model selection and retained
sequence/taxon identities. Mapping counts and improved paired coverage are not
yet claimed. Once the mapping completes, export its local PAE and extract/audit
native features before preparing paired sequence/3Di inputs. The existing
resource plan uses one CPU worker, 4 GB memory and 1 GB output headroom.

## Local residue mapping completed; native feature validation started

The local ESMFold mapping completed with all 5,121 models, 5,161 marker links,
201 taxa and 900,303 matrix-residue links. Seventy-four markers have at least
four linked taxa before confidence/coverage filtering. Readback checked every
mapping artifact hash, exact link equality with the independent prediction
audit, and every mapped amino acid against both the original NPZ sequence and
retained matrix character. Each mapped confidence agrees with NPZ confidence
within PDB decimal rounding; marker/taxon/column identities are unique.

Local PAE export has started under `results/structural_pae/esmfold-partial-v1`.
Pinned Foldseek extraction has finished under
`results/structural_alphabet/native-esmfold-partial-v1`, and coordinate-feature
reconstruction/validation is now running under
`results/structural_alphabet/coordinate-esmfold-partial-v1`. Extracted states
are not accepted for inference until that audit passes; joint PAE qualification
and paired-input preparation remain pending.

```bash
python scripts/export_local_marker_pae.py \
  --snapshot results/structural_markers/esmfold-partial-v1 \
  --output results/structural_pae/esmfold-partial-v1
python scripts/extract_structural_alphabet.py \
  --snapshot results/structural_markers/esmfold-partial-v1 \
  --output results/structural_alphabet/native-esmfold-partial-v1
OPENBLAS_NUM_THREADS=1 python scripts/audit_3di_features.py \
  --native results/structural_alphabet/native-esmfold-partial-v1 \
  --snapshot results/structural_markers/esmfold-partial-v1 \
  --output results/structural_alphabet/coordinate-esmfold-partial-v1
```

Mapping receipts/readback, native configuration and prelaunch resource estimates
are versioned as `metadata/esmfold_mapping_*`, `esmfold_native_extraction_*`
and `esmfold_coordinate_audit_resource_plan.json`. Existing alphaFold and local
snapshots remain separate, preserving prediction-source effects for assessment.

## Local confidence qualification and coverage contribution

The 5,121-model frozen local snapshot completed native coordinate validation
(1,183,341 residues; 1,173,099 valid structural states). All directional PAE
matrices exported losslessly. Joint six-residue pLDDT >=70 and maximum
context PAE <=10 retained 842,713 whole-model states. Independent NPZ readback
verified every artifact hash, model identity, PAE validity mask and count.

The same profile alignment and eligibility rule used for AlphaFold yielded
72 ready local marker alignments, 4,669 usable taxon-marker combinations and
199 represented taxa. All emitted AA/3Di identities, dimensions, missingness
masks and per-taxon coverage counts passed independent FASTA readback.
Local coverage includes 26 usable markers for Amoeboaphelidium occidentale
and 25 for Basidiobolus ranarum. The former adds an otherwise unrepresented
manifest lineage to the paired coverage inventory.

Comparing source-specific eligibility gives 458 entries with usable data from
at least one source (440 fungal entries and 18 outgroups), versus 322 in the
AlphaFold snapshot. Local predictions add 136 newly represented entries.
The 13,565 AlphaFold and 4,669 local usable taxon-marker combinations do not
overlap, consistent with targeting missing predictions. This is a coverage
union, not a pooled alignment or supported phylogenetic result. Predictor
calibration, uneven depth, copy resolution and uncertainty remain required.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/qualify_native_pae.py \
  --coordinates results/structural_alphabet/coordinate-esmfold-partial-v1 \
  --pae results/structural_pae/esmfold-partial-v1 \
  --output results/structural_alphabet/audited-esmfold-partial-v1
python scripts/prepare_paired_phylogenetic_inputs.py \
  --encodings results/structural_alphabet/audited-esmfold-partial-v1 \
  --snapshot results/structural_markers/esmfold-partial-v1 \
  --matrix results/phylogeny/profile-matrix-50-v1 \
  --output results/phylogeny/paired-inputs-esmfold-partial-v1
python scripts/summarize_paired_lineage_coverage.py \
  --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 \
  --output results/phylogeny/paired-lineage-coverage-esmfold-partial-v1
python scripts/compare_paired_source_coverage.py \
  --reference results/phylogeny/paired-inputs-gdm-expanded-v1 \
  --local results/phylogeny/paired-inputs-esmfold-partial-v1 \
  --output results/phylogeny/paired-source-coverage-v1
```

The source comparison also passed an identity check using the AlphaFold
snapshot on both sides: zero source-exclusive combinations, all 13,565
combinations shared, and 322 represented taxa. Summary receipts and tables
are versioned in metadata; raw structures and encodings remain outside Git.

## Same-sequence structural-alphabet predictor control

All 266 audited paired controls completed local mmCIF conversion, exact marker
mapping (270 links; 57,699 aligned residues), native Foldseek extraction,
coordinate reconstruction and directional PAE qualification. The comparison
binds the local native provenance to the exact completed control inference
configuration and requires identical full protein sequences to the selected
AlphaFold references. It uses existing qualified AlphaFold encodings.

`compare_predictor_alphabets.py` measures state disagreement at corresponding
residue positions under 12 confidence combinations (six-residue pLDDT 0/70/90;
context PAE unfiltered/5/10/15). Invalid native states are excluded. Descriptive
protein summaries require at least 50 residues and half the full protein in
both predictors. All 3,192 rows and the 400-cell state confusion matrix passed
count and threshold-nesting checks; three focused unit tests cover invalid
terminals, joint confidence, partner strata, sequence identity and coverage.

At pLDDT 70 and PAE 10, 187/266 controls meet coverage; 79 are excluded. Median
state disagreement is 0.14328 and median partner-change fraction is 0.18182.
For these same 187 proteins, unfiltered median state disagreement is 0.21481;
the median per-protein change after filtering is -0.07097. The retained residue
sets still change, so this is not a fixed-site effect estimate.

Within the filtered cohort, median state disagreement is 0.06587 for sites
with the same spatial partner and 0.48148 for sites whose partner changes.
These conditional descriptions do not establish that partner changes cause
disagreement. The proteins were deliberately selected across confidence,
length and lineage strata and are not a random estimate of the whole atlas.

Identical amino-acid sequences can therefore yield appreciably different
native structural states under these two predictors. This is a concrete
source-effect control, not evolutionary substitution, experimental accuracy,
or a universal noise floor. Pooling sources without sensitivity analyses could
confound branch interpretation. The ongoing local and AlphaFold branch fits
remain separate; leading results require matched-source and other controls.

Outputs: `results/prediction_controls/alphabet-comparisons-v1`; versioned
receipts and comparisons: `metadata/predictor_control_alphabet_*`.

```bash
python scripts/compare_predictor_alphabets.py \
  --inputs data/prediction_inputs/predictor-controls-v1 \
  --reference results/structural_alphabet/audited-gdm-expanded-v1 \
  --local results/structural_alphabet/audited-esmfold-controls-v1 \
  --output results/prediction_controls/alphabet-comparisons-v1
python -m unittest discover -s tests -p test_predictor_alphabets.py
```

Upstream control encoding stages reuse `convert_esmfold_snapshot.py`,
`map_marker_structures.py`, `export_local_marker_pae.py`,
`extract_structural_alphabet.py`, `audit_3di_features.py` and
`qualify_native_pae.py` with the corresponding `esmfold-controls-v1` paths.
Their completed provenance receipts and prelaunch resource plan are versioned.


### Longer marker queue: 513–768 residues

Prepared the complete canonical 513–768-residue band from the disjoint original
and follow-on marker candidate queues. It contains 5,512 sequences linked to
5,549 taxon/marker records across 279 taxa and 108 markers. A new frozen
exact-sequence reuse check found two sequences with verified available models,
leaving 5,510 prediction candidates. Complete protein sequences are retained.
The source inventory also records 7,883 canonical proteins above 768 residues
and 110 noncanonical proteins; this band does not resolve those gaps or the
full-proteome atlas.

Reproduce preparation with:

```bash
python scripts/prepare_long_marker_predictions.py \
  --inputs data/prediction_inputs/markers-v1 data/prediction_inputs/markers-followon-v1 \
  --predictions results/predictions/esmfold-marker-v1 results/predictions/esmfold-followon-v1 \
  --output data/prediction_inputs/markers-513-768-v1
```

The current run freezes 2,408 timing/memory receipts for proteins of 385–512
residues. Scaling runtime by length squared or cubed, using median and 90th
percentile normalized observed times, gives 41.82–66.70 inference GPU-hours.
The scheduling allowance is 100.05 hours (1.5 times the largest scenario).
Quadratic scaling of total allocated GPU memory gives 37.86 GB at 768 residues.
These are extrapolated scenarios, not measurements, complexity bounds or
confidence intervals. Allocator overhead, fragmentation and longer-sequence
behavior remain unmeasured. Plan one existing GPU, four CPU threads and 100 GB
output space. Predictions require a distinct configuration/output directory,
rechecked reuse and available GPU capacity, retaining the producer's OOM stop.
No longer-protein execution or controller was launched in this step.

Full readback verified every candidate sequence, the selected grid and all
5,549 source links within their originating candidate queue; available model
hashes; every resource observation; and all scenario calculations. An initial
readback expectation mistakenly counted all-marker records outside their
originating candidate queue; correcting that expectation removed duplicate
follow-on provenance without changing production inputs. Receipts and resource
plans are tracked in `metadata/long_marker_queue_receipt_v1.json`,
`metadata/long_marker_queue_readback_v1.json` and
`metadata/long_marker_resource_plan_v1.json`.


### Longer-marker execution queued after experimental controls

`scripts/advance_long_marker_predictions.py` is now running. It waits for the
exact experimental-control controller PID/start-time identity, requires its
successful receipt and a clean prediction chunk, then waits for the same GPU
to have no compute processes. Neither existing prediction job is interrupted.
Immediately before execution it reruns queue preparation into a new immutable
`data/prediction_inputs/markers-513-768-launch-v1` snapshot to refresh available
exact-sequence models and observed resource scenarios. It checks that the
selected length band, source receipts and complete marker links remain the
same; only the reusable/pending disposition may change.

The full remaining band will run in
`results/predictions/esmfold-markers-513-768-v1`, using the existing pinned
ESMFold checkpoint, batch size one, four CPU threads and maximum length 768.
The new output has its own prediction configuration. The controller requires
100 GB free storage and an unoccupied existing 49,140-MiB GPU; the extrapolated
memory estimate is not a guarantee of successful allocation. The existing
producer stops on OOM, and unsuccessful or interrupted completion is recorded
for review rather than automatically restarted. Final structure/PAE auditing
remains separate from clean execution accounting.

Controller configuration: `metadata/long_marker_prediction_controller_config.json`.
Log: `logs/long_marker_prediction_controller_v1.log`.
Controller outputs: `results/predictions/long-markers-controller-v1`.
The script compiled, verified its initial pins and entered its expected waiting
state. Input refresh and longer-protein prediction execution remain pending.


## Full-queue artifact audits and conversion queued

The original short-marker queue (10,522 eligible sequences, including its 128
cached models) and follow-on queue (4,252 eligible sequences) each have a queued
full-output handoff. `scripts/advance_prediction_snapshot.py` waits for the exact
producer PID and start time to exit, requires zero remaining eligible sequences,
no interruption or OOM, and the expected completed count. Pinned scripts, input
receipts and configuration must match before either stage runs. Existing audit
or conversion outputs prevent automatic reuse or overwrite.

The artifact auditor now additionally checks candidate hash identities and
uniqueness, eligible versus deferred sequence counts, per-model receipt identity,
length and exact artifact names, and reported low-confidence fractions and PAE
limits against NPZ arrays. Four snapshot tests (including corrupted metadata and
inconsistent remaining counts) and three link identity tests pass. Existing
coordinate, sequence, confidence, PAE and artifact-checksum checks remain active.

Controllers and configuration:

- `metadata/esmfold_original_snapshot_controller_config.json` →
  `results/predictions/original-snapshot-controller-v1`.
- `metadata/esmfold_followon_snapshot_controller_config.json` →
  `results/predictions/followon-snapshot-controller-v1`.

After successful full readback, all audited models are converted to sequence-
explicit mmCIF with atom-field roundtrip checks. Targets are respectively
`results/structures/esmfold-original-complete-v1` and
`results/structures/esmfold-followon-complete-v1`. Residue mapping, PAE binding,
native-feature qualification and evolutionary analyses remain subsequent steps.
Each handoff uses one CPU thread, an 8 GB memory planning allowance and 100 GB
output headroom; 0.2–12 hours is a conservative runtime planning envelope. It
requires no GPU inference or paid resources. These are queued jobs; no full-batch
audit or conversion completion is claimed at launch.


## Follow-on prediction queue completed; ecology queue launched

All 4,252 eligible follow-on sequences completed with zero cached models, zero
OOM events, no interruption and zero eligible sequences remaining. The producer
exited successfully. The 3,547 length/alphabet-deferred sequences remain explicit;
completion applies to the eligible short-protein queue, not the full atlas.
Execution records: `metadata/esmfold_followon_full_prediction_chunk.json` and
`metadata/esmfold_followon_full_prediction_config.json`.

The full artifact auditor was observed live reading all 4,252 model receipts.
Audit results and subsequent mmCIF conversion remain pending. Independently,
the ecology controller verified the same completed chunk and launched all 675
queued ecology sequences on the freed GPU. Its launch receipt is
`metadata/esmfold_ecology_prediction_launch.json`; the inference process was
observed live. Neither launch establishes new confidence-qualified marker
coverage or an ecological association result.


The full follow-on artifact audit passed for 4,252 models, 4,304 marker links,
84 taxa and 98 markers. Independent exported-table readback matched the entire
eligible FASTA identity set and all originating marker-link rows. Median length
is 309 residues; median model mean CA pLDDT is 85.878, which is not experimental
accuracy. The full input queue spans more taxa than its short-protein eligible
subset. Conversion is now running; residue mapping and confidence-qualified
paired coverage remain pending. Receipt:
`metadata/esmfold_followon_full_audit_receipt.json`.
