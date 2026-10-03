# Complete independent native alignment replay

The [35-pin full plan](../metadata/independent_native_alignment_replay_plan_20261002.json)
launched October 2 at 23:13 EDT. It retains all 1,620 selected original chain
identities and all 405 groups: 1,618 intact chains, 403 complete quartets and
two failed chains in two unresolved quartets. All six intact chains in those
unresolved quartets are included. The producer is running; its complete raw
serialized readback and provenance closure are queued. This does not qualify
an ancestral posterior or a new ancestral structure.

## Exact scope and source lineage

The [completed source inventory](../metadata/independent_native_alignment_inventory_20261002.json)
and its [original successful journal evidence](../metadata/independent_native_alignment_inventory_completed_20261002.json)
fix the selected complete attempts, chain IDs, seeds, input alignments, sample
audits, runtime trees, candidate nodes, stored state arrays and coordinate
manifests. No original attempted samples are concatenated or overwritten.
Both native failures retain their original dispositions.

Full replay covers 163,418 saved native alignments, 653,672 candidate-node
frames, 9,732,673,504 projected state observations, 96,363,104 per-chain
anchor coordinates and 708,466 unanchored-residue observations. Cutoffs 250
and 500 each contribute 2,119,988,288 count-array cells, totaling
4,239,976,576 cells across all 22 amino-acid/unknown-X/gap states. Repeated
anchors are not independent biological sites or evolutionary observations.

The [source loader](../scripts/independent_native_alignment_replay_sources.py)
requires the closed full recovery/diagnostic graph and the independently
decoded inventory. Every inventory entry is reconciled with the selected
native attempt and all original state/coordinate receipt links. Full source
hashes are checked before and after each stage. The inventory contains
50,066 bindings, totaling 65,141,744,592 bytes per full hash pass; added new
stage bindings are also checked. The source inventory's input FASTA/tree
decoding alone was not raw saved-sample verification.

## Native decoding and exact array comparisons

The [manual implementation](../scripts/independent_native_ancestral_alignment.py)
imports neither Biopython nor the original saved-block parser, anchor builder
or projection. For every chain, the [replay engine](../scripts/independent_native_alignment_replay.py)
decodes the input FASTA and runtime Newick labels/tips and checks the complete
coordinate schema, listed tip lengths, node order, input SHA and declared
coordinate ordering. Each native stream must contain exactly iterations
0–1000 in increments of ten, with no missing, duplicate or reordered samples.
Every saved sample must contain the exact runtime node set and equal-width
valid sequences. Concrete extant residues cannot change; expected input X
residues retain wildcard semantics.

Independent projection selects ancestral letters at every sorted tip's
nongap columns, retaining duplicated anchors. Candidate residues in columns
with no extant residue contribute to unanchored counts. Every candidate
ungapped length must match the original per-draw audit. All projected state
bytes and all unanchored counts must match the original arrays exactly.

Freshly decoded complete trajectories supply both cutoff count arrays.
Histograms use integer coordinate/state keys and `bincount`, in blocks of
4,096 coordinates, separately from the original per-label comparisons.
Every count cell, shape, dtype and retained-draw count is checked. Both
cutoffs retain the original strict greater-than iteration rule: 75 and 50
saved samples per chain. There is no numerical tolerance for these discrete
comparisons and no trajectory thinning or state omission.

Each logical sample is serialized with its iteration, width, full alignment
digest, projected-state digest, candidate lengths, unanchored counts and an
explicit false scientific-eligibility flag. Whole-chain checkpoints include
full-array/count digests, every original source identity and count totals.
Source-to-runtime clade identity inherits the closed audit; the manual tree
decoder checks grammar and labels, not a new biological topology proof.

## Checkpoints, complete readback and software evidence

The [producer](../scripts/prepare_independent_native_alignment_replay.py)
streams logical frame JSONL, flushes/fsyncs and atomically publishes it only
after every sample and count array passes. Immutable whole-chain checkpoints
cover intact and failed chains. Interruption recovery re-decodes all saved
frames and compares saved rows without rewriting accepted files. Completed
producers refuse restart, and a stage lock prevents duplicate writers.

The [reader](../scripts/readback_independent_native_alignment_replay.py)
re-decodes every original native alignment, rebuilds every projection/count
array and compares all saved frame rows and checkpoints. It rejects missing,
extra, changed and foreign artifacts, reconciles every full-scope count and
rechecks source/artifact hashes. It shares the new decoder and is not a third
independent implementation. Final accounting closure requires complete full
artifact/source proofs and both exact original producer/reader completion
journals. A live handle or default collected-service success field alone
does not prove completion.

The [full workflow tests](../metadata/independent_native_alignment_replay_grid_validation_20261002.json)
retain all 1,620 synthetic chains/405 groups and two explicit failures, including
all six intact chains in unresolved quartets. Four variable alignment templates
include lowercase/multiline FASTA, reordered records, tip X wildcards, shifting
gap positions, all 22 candidate states and unanchored residues. Expected arrays
come from the original projection, not the new decoder. The complete software
grid replays 163,418 saved samples/653,672 candidate frames, 2,614,688 state
observations and 1,139,072 count cells. Six histogram fixtures cover both
retained-draw counts and coordinate-block boundaries.

Every serialized row was reconstructed; interrupted checkpoints stayed
byte-identical, and completed producer/alternate-reader restarts were refused.
Fourteen rehashed false exports were rejected: missing/duplicate frames,
changed alignment/projection digests, counts, candidate lengths, scientific
claims, cutoff/full-state digests, seeds, false failure dispositions, hidden
cutoffs and omitted failed/intact-unresolved chains. Eight original-source
alterations were also rejected: changed states/unanchored arrays, either
cutoff histogram, saved iteration arrays, duplicate native iterations,
changed observed residues and missing runtime nodes. These are software
contracts with synthetic source journals, not a biological pilot or complete
production verification.

## Execution, environment and resources

The [launch inventory](../metadata/independent_native_alignment_replay_launches_20261002.json)
records original producer controller 2716502, reader 2716506 and closure 2716513.
The [first execution observation](../metadata/independent_native_alignment_replay_execution_checkpoint_20261002.json)
verified exact live original identities, the native producer and installed
limits. No chain checkpoint was complete at that observation; the producer
was checking sources. These are original new invocations, not restarts of the
native sampler or failed categorical run.

The [complete resource inventory](../metadata/independent_native_alignment_replay_resources_20261002.json)
records 45,530,706,578 native alignment bytes. The largest native file is
1,392,681,456 bytes and the largest original state array is 115,136,768 bytes.
One chain is resident at a time. The estimated array workspace is
1,058,881,024 bytes; it is not a hard bound. Installed per-stage cgroups enforce
two CPU equivalents, 32 GiB memory and no swap, with one BLAS thread. Output
planning allows 4 GiB and requires 104 GiB free disk. The 1–48-hour range per
stage is uncalibrated planning, not a project or convergence ETA. No GPU,
new native inference or paid resources were launched.

The [environment definition](../environments/independent-native-alignment-verification-20261002.yml)
records Python 3.10.13, NumPy 2.2.6, SciPy 1.15.3, psutil 7.2.2 and
Biopython 1.86. Biopython is used by the original expected-state fixture
module, not the new production decoder. Reproduce tests and inspect the
existing run without relaunching production:

```bash
python scripts/check_independent_native_alignment_replay.py --output NEW_GRID_CHECK.json
python scripts/record_independent_native_alignment_replay_checkpoint.py \
  --output NEW_EXECUTION_OBSERVATION.json
```

Full production replay/readback/provenance closure, independent clade mapping,
adequate joint sampling, biological model/root/topology/predictor qualification
and all eight evolutionary aims remain incomplete.

At 23:18 EDT, a [second actual-handle/cgroup observation](../metadata/independent_native_alignment_replay_execution_checkpoint_20261002_v2.json)
confirmed 20 complete producer chain checkpoints, 2,020 native alignments,
502,223,712 projected states and 218,790,528 cutoff-count cells, with 8,946
unanchored observations. Full source checks passed before these decodes.
These are unclosed producer progress; complete serialized and provenance
verification remains pending.


## October 3 producer completion

The original producer terminated successfully at 00:14 EDT after the full
1,620-identity scope and final source checks: 1,618 intact chains, both
failures, 163,418 alignments, 9,732,673,504 states and 4,239,976,576
cutoff-count cells. The [latest original-handle observation](../metadata/independent_native_alignment_replay_execution_checkpoint_20261003_v4.json)
verified its actual successful terminal journal; the separate original
serialized reader remained live, reaching 872/1,620 chains by 00:35 EDT.
Final serialized/source/artifact/journal closure remains pending. No
ancestral ensemble was qualified. The separate node-correspondence stage
has closed; a [full sampling resource/numerical census](baliphy-horizon-resources-20261003.md)
now informs the longer-horizon design without launching new native work.


## October 3 complete serialized and provenance closure

The [full completion](../metadata/independent_native_alignment_replay_completed_20261002.json)
now verifies all 163,418 saved alignments, 653,672 candidate frames,
9,732,673,504 anchored state observations, 708,466 unanchored observations
and 4,239,976,576 cutoff-count cells. All 1,618 intact chains and both
failed identities remain included, including the six intact members of
unresolved quartets. Complete serialized reconstruction and provenance
closure bound 53,331 source/artifact hashes and both original journals.

The [final actual-handle observation](../metadata/independent_native_alignment_replay_execution_checkpoint_20261003_final.json)
at 05:27 UTC rechecked the compact archive hash and verified successful
terminal journals for the original producer, reader and closure controllers
2716502/2716506/2716513. This observer did not repeat the whole large-data
hash/decode stage; that work is recorded in the completed archive.
The reader shares the new decoder and is not a third independent algorithm.

Native source/serialization integrity and the separately closed clade mapping
are established for the complete selected scope. They do not qualify posterior
mixing, the supplied biological root, model adequacy, predictors or any of
the eight biological aims. The exact categorical ESS discrepancy remains
numerical review, and both earlier allocation failures remain unresolved.
The [new full reference-startup preflight](baliphy-reference-initialization-20261003.md)
is separate and does not repeat or concatenate these native samples.
