# Joint sequence/structural directions in original gene-tree contexts

The full integration's producer has completed **283,409 original gene-tree
contexts, 566,818 design records, 214,461 tied references and 428,922 logical
sides**. It links all 27,056 source-ready ordered physical triples and 730,512
joint direction groups back to their original parent and reference identities.
Full independent readback and closure remain pending. This extends the
[completed structural-context directions](full-triad-context-contrasts-20261002.md)
to the [joint sequence/structural physical control](full-triad-joint-directions-20261002.md).

## Full grid and count units

All 27 full/pLDDT70/both-mask × FAMSA/MAFFT/both-method ×
reference-common/cycle-consistent/both-core scenarios and all six original
quality screens remain explicit. The physical source retains separate and
deduplicated joint envelopes, complete order grids, unavailable/nonunique
states and numerical uncertainty, without favorable correspondence selection.

| Full integration output unit | Count |
| --- | ---: |
| Original contexts | 283,409 |
| Availability/sequence-first design records | 566,818 |
| Tied references across those designs | 214,461 |
| Logical duplicate/reference sides | 428,922 |
| Reference × scenario × screen cells | 34,742,682 |
| Context × design × scenario × screen cells | 91,824,516 |
| Repeated five-policy direction decisions | 459,122,580 |
| Disjoint ten-state count rows | 32,400 |

These are dependent descriptions of overlapping contexts and settings, not
independent evolutionary events or multiple-testing hypotheses.

## Original source gates and fixed references

Each record embeds the unchanged original work design: guide, family, taxon,
gene-tree node, duplicate genes, model IDs/versions, all reference genes and
their original lexical/tied order. Both availability and sequence-first designs
are retained. A missing lexical model is never replaced by a later modeled tie.

Every reference has three separately reconstructed source gates: readiness,
native membership in its guide, and native membership in both guides. Readiness
requires an eligible parent, distinct versioned models and distinct model IDs,
the original duplicate disposition, and all three pairs in the unchanged
measurement designs. A physically measured triple cannot override a rejected
parent or a same-ID/different-version context.

References link to normalized `[triad_id, mask, sequence_method, core]` keys for
all 27 scenarios. Scenario identifiers explicitly join those three axis labels
with `|`; the ordered `scenario_axes` list in the plan and receipt defines their
mapping. The output uses this flat scenario dimension rather than hiding method
selection inside a core label. Each scenario retains six screen cells and all
three source-gate direction states. Missing or gated direction is null; it is
distinct from a measured zero/near-zero contrast. A physical group can be
present but fail every quality screen, including missing or nonunique fit cases.

## Reference pools and direction

Five fixed policies retain the original lexical reference under each of the
three source gates, all eligible native-both ties, and all original native-both
ties. The any-eligible policy compares **every eligible tied reference** and
retains three gate-specific eligibility counts plus four native-both direction
support counts. Agreement among an incomplete eligible subset does not become
complete original-tie agreement. The all-original policy requires a nonempty
tie set and every original reference eligible.

Quality eligibility and direction consistency remain separate. Differing
qualified reference categories produce `reference_direction_disagreement`,
which can include positive versus sign uncertain and need not mean opposite
signs. Positive, negative, within numerical tolerance, sign uncertain, reference
disagreement, no ties, source exclusion, geometry exclusion, no eligible reference
and incomplete all-tie eligibility form ten disjoint summary states. No reference,
method, mask, core or order is reselected to obtain a desired direction.

The inherited 1e-9 Å threshold is a numerical classification boundary, not an
effect-size cutoff or a confidence interval. Sequence-derived correspondences
use the same predicted coordinates and do not remove predictor circularity.

## Independent verification, restart and resources

The producer links closed joint-direction groups to the original contexts.
The reader independently reconstructs raw original model-role SHA identifiers,
distinct-model/pair/parent/native gates and every scenario key. SQLite separately
reconstructs lexical directions, all eligible pools, all-tie completeness,
direction support, counts and source denominators. Every exported record, scalar
type, source-row order and deterministic checkpoint must agree. Complete source/
artifact hashes and both original invocation-linked completion/resource journals
gate final closure. The reader does not import producer projection/policy logic.

Software contracts passed 30 contexts, 64 reference ties, eight physical triples,
216 joint physical groups, all 27 scenarios, 8,100 policy decisions and 5,400
one-screen summary rows. Cases cover rejected measured parents, missing lexical
models, empty pools, model-version exclusions, native-guide disagreement,
reference-category disagreement, incomplete eligible pools, mask/core/method
direction changes, unavailable/nonunique fits and near-zero measurements.
All 23 rehashed false exports were rejected. Interrupted recovery rechecked one
committed four-context chunk; completed restart was refused. Production source-
closure I/O is stubbed only in software fixtures. The first fixture attempt
failed to write its receipt because its output directory did not exist; after
adding directory creation, the full suite passed at a new v2 path.

An exclusive lock is held through atomic receipt creation. Atomic configuration
and deterministic 1,000-context chunks produce 284 full-design checkpoints.
Checked interrupted recovery regenerates committed chunk bytes and rebuilds
exports; completed outputs remain immutable. Reproduction requires new plan/
output/unit identities or a verified terminal original process before restart.

Prelaunch resources: two CPU cores, 32 GiB RAM, no swap, one BLAS thread,
192 GiB output allowance and 100 GiB free-disk reserve. The measured completed
nine-scenario contextual stage occupies 1,554,399,112 bytes including its full
hash archives; this stage triples the scenario grid and indexes 730,512 physical
groups. The conservative allowance also covers uncompressed export accounting
and duplicate compressed checkpoints. The 1–24-hour producer/reader planning
range is uncalibrated and excludes dependency waits; it is not a completion ETA.
No GPU, paid resources or changes to existing scientific jobs are involved.

- [Frozen full plan](../metadata/full_triad_context_joint_directions_plan_20261002.json)
- [Prelaunch resource estimates](../metadata/full_triad_context_joint_directions_resources_20261002.json)
- [Passed software contracts](../metadata/full_triad_context_joint_directions_fixture_validation_20261002.json)
- [Original launch identities](../metadata/full_triad_context_joint_directions_launches_20261002.json)
- [Closed-source loader](../scripts/full_triad_context_joint_direction_sources.py)
- [Producer](../scripts/project_full_triad_context_joint_directions.py)
- [Independent raw-gate/SQL reader](../scripts/readback_full_triad_context_joint_directions.py)
- [Software checks](../scripts/check_full_triad_context_joint_directions.py)
- [Launcher](../scripts/launch_full_triad_context_joint_directions.py)

The original producer, reader and closure wait for the original joint-direction
closure, whose launch identity is pinned by the plan. Large outputs stay outside
Git in `results/structural_comparisons/full-triad-context-joint-directions-20261002-v1/`.
The future completion locator is
`metadata/full_triad_context_joint_directions_completed_20261002.json`.

Reproduce the software contracts with a new output path:

```bash
python scripts/check_full_triad_context_joint_directions.py --output results/software-checks/context-joint-directions-new-run/receipt.json
```

For fresh production, prepare new pinned plan/output/launch identities after
required source closure. The entry points are:

```bash
python scripts/project_full_triad_context_joint_directions.py --plan metadata/new_context_joint_direction_plan.json
python scripts/readback_full_triad_context_joint_directions.py --plan metadata/new_context_joint_direction_plan.json --output results/structural_comparisons/new-context-joint-directions/readback.json
```

Run the reader after successful production. Existing live jobs and completed
outputs must not be duplicated or overwritten.

This is a reference/correspondence control supporting duplication-associated
comparisons. It does not establish accepted gene/species phylogenies, ancestral
polarity, structural acceleration, selection or calibrated duplication effects.
Actual matched backgrounds, predictor/domain/PAE controls, accepted species and
reconciled gene framework, family/taxon dependence, missingness/sampling and
statistical calibration remain required. All eight scientific aims remain
incomplete; GPU prediction remains paused.
## Current production checkpoint — October 2 UTC

The original producer completed all 283,409 contexts/214,461 ties and
459,122,580 dependent policy decisions. Its original process is absent and its
invocation-linked completion/resource journal passed the [fresh full runtime
check](../metadata/project_runtime_checkpoint_20261002_v20.json). Independent
original-source/SQL readback and final closure remain live; there is no final
completion claim for this context stage. The upstream [joint physical directions](full-triad-joint-directions-20261002.md#complete-production-and-published-results-october-2-utc)
are now fully closed and published. This updates the earlier queued status.
All original gates remain intact, all eight aims incomplete, GPU prediction paused.
