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
