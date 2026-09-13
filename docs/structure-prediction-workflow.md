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
