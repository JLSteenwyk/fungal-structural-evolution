# Independent source-to-runtime ancestral node correspondence

The [49-pin complete plan](../metadata/independent_native_clade_mapping_plan_20261002.json)
launched October 2 evening. It retains all 1,620 original chain identities,
405 model groups and 135 effective input groups. All 1,618 intact chains,
including six in unresolved quartets, are checked; both failed chains retain
their original dispositions. Full production, serialized reconstruction and
[provenance closure](../metadata/independent_native_clade_mapping_completed_20261002.json)
are complete. The supplied biological root, tree and
model remain conditional assumptions.

This extends the [full native residue replay](independent-native-alignment-replay-20261002.md).
That replay checks every saved state/count but inherits the original
source-to-runtime candidate mapping audit. This separate stage independently
reconstructs the topology and candidate correspondences from the original
source and same-attempt runtime trees. No running native-residue or categorical
stage was edited or restarted.

## Tree, branch and candidate checks

The [manual tree implementation](../scripts/independent_native_ancestral_topology.py)
uses the separately validated lexical tokenizer, an iterative named-Newick
grammar and postorder integer descendant bitsets. It imports neither
Biopython nor the original descendant-index function. Source and runtime
tips must agree exactly. Every named node has a unique descendant set;
duplicate names, repeated leaves, overlapping child sets and ambiguous
unary clades are rejected. All source and runtime rooted clade sets must
match, including the full-tip root. Child ordering and internal runtime
labels can change without changing clade identity.

Every nonroot edge length is parsed as a finite decimal and compared by
exact rational subtraction. The original strict absolute threshold remains
`difference < 1e-10`. Missing nonroot lengths and differences at or above
the threshold fail. Root edge lengths are excluded, as in the original
audit; their serialized values remain available. Negative branch pairs
retain a separate review count. This is correspondence of conditional
phylogenetic branch parameters, not dates or physical structural displacement.

The [full correspondence engine](../scripts/independent_native_clade_mapping.py)
reconciles all original IDs, seeds, selected attempts, input trees/alignments,
runtime trees and sample-audit links with the closed full inventory. It
independently selects the original profile-guide whole/domain mapping rows
for each family, reconstructs retained-tip sets, checks their source-node
clades, derives the runtime labels and compares the saved mapping. Levels,
labels, iteration identities and uniqueness are checked for all 404 candidate
audit entries per intact chain. All four candidates remain distinct, with
level 3 spanning the conditional full-tip root.

Production scope is 161,156 source/runtime node pairs, 159,538 nonroot branch
pairs, 6,472 candidate-node pairs and 653,672 saved-audit candidate mappings.
All 1,618 assumed-root candidates retain `biological_root_accepted: false`.
This stage does not independently decode the saved alignment sequences or
recompute their residue/count arrays; those are handled by the separate
native replay. Correct node correspondence does not prove a biological root,
model adequacy, joint mixing or a qualified ancestral ensemble.

## Serialization, complete tests and provenance

The [producer](../scripts/prepare_independent_native_clade_mapping.py)
publishes immutable whole-chain checkpoints with every node/clade/length pair,
ordered tips, candidate masks/tip digests and original source identities.
All descendant masks use the checkpoint's sorted full-tip order. Interrupted
checkpoints are rebuilt from raw source/runtime trees and compared without
rewriting accepted files. Completed producers refuse restart.

The [serialized reader](../scripts/readback_independent_native_clade_mapping.py)
reconstructs every original tree, node pair, length and candidate identity,
compares complete checkpoints and manifests, reconciles all scope partitions
and rechecks full source/artifact hashes. It shares the new tree decoder and
is not a third independent tree implementation. Final closure requires the
complete proof graph and both exact original producer/reader completion
journals. Missing closure is not inferred to be a job failure.

The [complete software workflow](../metadata/independent_native_clade_mapping_grid_validation_20261002.json)
passed all 1,620 synthetic chain entries/405 groups/135 inputs, including
both failures and six intact-unresolved chains. It checks 14,562 synthetic
node pairs, 6,472 candidate pairs and 653,672 candidate-frame mappings.
Forty randomized rooted tree correspondences matched the original Biopython
descendant index. Quoted labels/comments, a 2,399-node iterative tree,
negative-branch review and both sides of the strict branch threshold were
checked. Eight invalid correspondences were rejected.

Interrupted checkpoints stayed byte-identical; completed producer and
alternate-reader restarts were refused. Twelve rehashed false exports were
rejected: missing/duplicate nodes, changed clades/lengths/candidate labels,
false root acceptance, foreign seeds, wrong frame counts, hidden failures
and omitted failed/intact-unresolved chains. Four altered source audits were
also rejected: wrong candidate levels/runtime labels and duplicate/missing
frame mappings. The repeated source fixtures and their journal fields are
synthetic, not a biological pilot or actual production provenance.

## Original execution and resources

The [launch inventory](../metadata/independent_native_clade_mapping_launches_20261002.json)
records producer controller 2735141, reader 2735145 and closure 2735152.
Both prerequisite original inventory/diagnostic completion journals were
checked before launch. The [first execution observation](../metadata/independent_native_clade_mapping_execution_checkpoint_20261002.json)
verified all three original live PID/create/command identities, the native
producer and installed two-CPU/32-GiB/no-swap limits. No checkpoint was yet
complete while initial source hashes were being checked.

The [complete resource inventory](../metadata/independent_native_clade_mapping_resources_20261002.json)
covers 90 source-tree paths/293,355 bytes and 1,618 runtime trees/3,384,640
bytes, with at most 1,243 runtime nodes per chain. The full inherited
source graph has 50,066 bindings/65,141,744,592 bytes per hash pass, and
large-data hashes are checked afresh before and after each stage. New
stage bindings are also included. The two-GiB object workspace is an
estimate, not a bound. Output planning allows four GiB and requires 104 GiB
free disk. The 0.1–12-hour range per stage is uncalibrated planning, not a
finish or convergence ETA. One BLAS thread is used. No native inference,
GPU prediction or paid resources were launched.

Use the [native verification environment](../environments/independent-native-alignment-verification-20261002.yml).
Reproduce the software checks or inspect the existing run:

```bash
python scripts/check_independent_native_clade_mapping.py --output NEW_GRID_CHECK.json
python scripts/record_independent_native_clade_mapping_checkpoint.py \
  --output NEW_EXECUTION_OBSERVATION.json
```

The [October 3 execution check](../metadata/independent_native_clade_mapping_execution_checkpoint_20261003.json)
verified all three original successful terminal journals and the completion
archive hash. Full closure binds 51,725 source/artifact hashes and the two
original producer/reader journals. All 161,156 node pairs, 159,538 nonroot
branch pairs and 653,672 candidate mappings passed. Both failed-chain
dispositions remain unchanged; no negative branch pairs were observed.

Do not relaunch the three existing production invocations. Full
residue/categorical closure, adequate joint
sampling, accepted root/tree/model/predictor controls and all eight biological
aims remain incomplete.
