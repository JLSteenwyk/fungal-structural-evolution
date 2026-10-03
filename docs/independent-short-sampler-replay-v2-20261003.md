# Observation-aware full short-sampler replay, version 2

The corrected full replay is queued for **all 1,620 original roles, 405
quartets, 135 effective inputs and 324 configuration aliases**, behind the
original sampler's complete source/artifact/journal closure. It checks both
separate tip draws against the observed input, retains discrepancies at
unknown residues, and never treats separate FASTA/category logs as one joint
trajectory. Every failed or invalid original disposition remains required.
No new native inference, GPU work or charges are launched.

## Required observation constraints

The [full input census](../metadata/independent_short_sampler_full_input_ambiguity_census_20261003.json)
found three `X` residues in one of 39 distinct alignment files, affecting 24
roles. A known observed amino acid must match both the saved FASTA letter
and the residue decoded from the category/state record. Both logs agreeing
on a different amino acid is still rejected. Array lengths, category/state
support, native alphabet order, rate matrix dimensions/values/normalization,
condition schema, iteration, tip labels and strict JSON rules remain checked.

An input `X` represents an unknown residue. Both native logs must contain
valid amino acids at that position, but their conditional draws may differ.
Version 2 retains each difference with the tip label, zero-based ungapped
source-tip offset and both letters. It does not overwrite either draw or
infer the category allocation for the separate FASTA draw. It marks
`tip_logs_joint_trajectory_available=false` even when the letters coincide.
Internal category labels remain unavailable in this historical horizon.

The residue projection arrays are reconstructed from the legacy FASTA only;
category/state records are checked and retained separately. Source/runtime
clades, rooted candidate correspondence, all three saved alignments and all
four candidate projections retain the [original replay scope](independent-short-sampler-replay-20261003.md).
Complete producer/reader value, dtype, axis, hash, export and artifact-set
comparisons remain required. No successful projection array is manufactured
for an unresolved original role, and intact roles in an unresolved quartet
are not promoted into an accepted quartet.

## Software evidence and its limits

The [software qualification](../metadata/independent_short_sampler_replay_v2_software_validation_20261003_v1.json)
passed:

- Complete producer/reader serialization with all 1,620 real role metadata
  records and explicitly mocked source admission/native results. This
  includes every one of the 24 ambiguity-affected roles, 216 repeated unknown
  positions and 72 artificial differing draws. Two artificial failed roles
  and all six intact members of their unresolved quartets remain explicit.
- Read-only decoding of three actual available native roles, one per prior:
  all nine alignments and 36 candidate frames match the original independent
  anchor-projection oracle.
- Read-only checks of nine retained native synthetic ambiguity frames across
  all three priors. Five unequal draws reproduce the old reader's rejection;
  version 2 accepts their observed constraints and retains both letters.
  These are actual saved software-test outputs, not new fungal inference.
- All 400 possible pairs of amino-acid draws at one unknown position, with
  disagreement counts retained for the 380 unequal combinations.
- Fifteen malformed property/JSON cases, three extra observation cases and
  fourteen restart, omission, serialization and false-acceptance alterations
  rejected. These include hidden disagreements, changed ambiguity counts and
  invented joint-trajectory claims.
- One capped pure native alphabet introspection probe, without MCMC or
  alignment sampling. The 22-state gap/X projection fixture still matches
  the separate original projection implementation.

The successful original software invocation was verified by its actual
process identity, invocation-linked completion/resource journal and 98 source
bindings. It used 11.971485 CPU seconds. The launch configured two CPUs,
16 GiB RAM and no swap, but the subsequent live cgroup capture occurred after
the fast unit was collected. Missing live limits and memory-peak observations
are explicitly unavailable; collected default fields are not substituted.
No software or native attempt was rerun to replace the observation.
[Actual completion evidence](../metadata/independent_short_sampler_replay_v2_software_execution_20261003.json).

Full-grid serialized software checks use mocks. The actual existing native
checks are limited read-only QC. Neither proves full production replay,
posterior mixing, independent likelihood correctness or biological adequacy.
The separate [future joint logger](baliphy-joint-node-logger-20261003.md)
has its own software evidence and prepared source grid; it cannot restore
unlogged historical ancestral category trajectories.

## Frozen production queue and resources

The [immutable version-2 plan](../metadata/independent_short_sampler_replay_v2_plan_20261003.json)
binds 1,008 sources. Its new output root is
`results/ancestral/full-independent-short-sampler-replay-20261003-v2`.
The original sampler, version-1 reader, native attempts, limits and schedules
remain untouched. Version 2 depends directly on verified sampler closure,
so a version-1 observation mismatch does not invalidate this new dependency.

The [launch inventory](../metadata/independent_short_sampler_replay_v2_launches_20261003.json)
captures original producer/reader/closer controllers 2915414/2915420/2915424.
Each has two CPUs, 16 GiB RAM, no swap and one-thread BLAS settings. Reader
and closer wait for their exact original dependencies and actual completion
journals. Completed stages refuse restart; partial or foreign artifacts must
be reviewed, not silently adopted. Full native/source/artifact hashes and
both original replay stage journals are required for final accounting closure.

Planning allows 8 GiB output, requires 64 GiB free disk and retains the
uncalibrated 0.05–8 active hours per stage. If all original roles pass, the
potential scope is 4,860 saved alignments, 19,440 candidate frames and
293,964,912 projected uint8 bytes; the largest single-role projection contains
3,419,904 bytes. These counts do not bound parser, native FASTA/JSON or
transient proof workspace. At preparation, 307.73 GiB RAM and 10,372.25 GiB
disk were available. This is a resource estimate, not a completion ETA.

At 08:45 UTC, all six original version-2 replay/sampler handles and live caps
were reverified. The sampler had 869 successful unclosed role checks;
version-2 replay had zero role checks and remained queued. The
[version-1 handles and 1,006 pins](../metadata/independent_short_sampler_replay_execution_checkpoint_20261003_v7.json)
were independently reverified. Both versions' production disposition/readback/
closure remain pending.
[Version-2 runtime evidence](../metadata/independent_short_sampler_replay_v2_execution_checkpoint_20261003_v1.json).
The resource observer remains live with
[fresh evidence](../metadata/baliphy_sampler_resource_observation_execution_checkpoint_20261003_v6.json),
and the broader covariance/timing queue remains
[verified](../metadata/full_shared_entity_timing_execution_checkpoint_20261003_0847.json).

Original allocation failures, adequate longer ensembles, root/model/predictor
controls, phylogeny/reconciliation/dating, calibration and all eight biological
aims remain open. GPU structure prediction remains paused.

## Reproduction and interpretation

```bash
python scripts/record_independent_short_sampler_replay_v2_checkpoint.py \
  --output NEW_V2_RUNTIME_CHECKPOINT.json
```

Observe the original handles; do not relaunch frozen stages. Software checks
may use fresh audit/receipt destinations with
`scripts/check_independent_short_sampler_replay_v2.py`; the native probe has
an 8 GiB address-space cap, 120 CPU seconds, 150 wall seconds and 1 MiB
per-file cap. The [software resource estimate](../metadata/independent_short_sampler_replay_v2_software_resource_plan_20261003.json)
and [existing Python environment](../environments/independent-short-sampler-replay-20261003.yml)
document those settings. Source/native inputs and executables remain hash
bound; large arrays, fixtures and evidence archives stay outside Git.

The [version-2 data dictionary](../metadata/independent_short_sampler_replay_v2_data_dictionary_20261003.tsv)
adds observed-constraint and disagreement fields to the
[original output dictionary](../metadata/independent_short_sampler_replay_data_dictionary_20261003.tsv).
The 20-state native category alphabet is `ARNDCQEGHILKMFPSTWYV`;
22-state projection arrays use `ACDEFGHIKLMNPQRSTVWYX-`. Repeated frames,
chain roles and tip anchors are dependent measurements, not independent
biological replicates.
