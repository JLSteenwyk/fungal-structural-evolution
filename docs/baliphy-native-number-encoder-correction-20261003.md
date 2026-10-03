# Corrected ancestral property encoding and full historical replay

The installed Haskell double formatter still truncates scientific exponents.
The new V5 logger bypasses it for substitution-model properties using the
installed native CJSON encoder. No installed binary or library was patched.
All failed runs, historical numbers and original model programs remain intact.
[Original diagnosis](baliphy-native-number-formatting-20261003.md).

## Exact scope of the correction

`baliphy_joint_node_logger_v5.py` composes the preserved full-node/same-record
sequence logger with two rendering edits: a qualified encoding import and
native CJSON encoding of the `properties` field. Reversing all five edits
recovers each original program exactly. The probability model, distributions,
priors, initializer, observation and ancestral draw expressions stay unchanged.
The future frame readers still require finite nonnegative rates with mean one
at relative/absolute tolerance `1e-10`; neither decoder was relaxed.

Full static qualification covers every **405 programs and 1,620 roles**.
Nine bounded native software runs cover all three priors. Three corrected
runs match retained same-seed scalar logs, runtime trees, legacy alignments
and nonproperty joint-frame fields byte for byte. Three small-alpha fixture
runs and three corresponding instrumented counterparts reproduce six old
normalization failures while the corrected frames pass. The test-only oracle
roundtrips all 80 native rate cells exactly in nine frames: **720 exact native
double roundtrips**. Adding the oracle does not change any saved sampler output.
Fixed alpha and oracle logging occur only in software fixtures.

Two independently implemented sequence/category readers validate all 18
corrected frames and 72 ancestral records. All candidate arrays undergo
serialized key, shape, dtype and value readback. Both readers reject all
108 malformed cases, including incorrect rate means. These are software
checks within the full design, not a biological pilot or adequate posterior.
[Native/static qualification](../metadata/baliphy_joint_node_logger_software_validation_20261003_v5.json).

Prepared future sources retain all 405 models, 135 effective inputs, 324
configuration aliases and four roles per model. Their 1,620 deterministic
seeds are disjoint from the 4,866 explicitly forbidden original, earlier
future and fixture seeds. Future models contain neither fixed-alpha changes
nor the oracle. [Source/role manifest](../metadata/baliphy_joint_node_logger_future_models_20261003_v5.json).

## Historical replay: sequence accounting with unqualified rates

The V3 historical reader retains the original logged rate cells unchanged.
It records the same strict normalization check as a quality flag and keeps
**all historical rates unqualified**, even if their mean passes. It does not
renormalize, infer missing exponents or claim that twelve failures are the
complete corruption census. Corrected future readers remain strict.

Full replay and separate serialized readback have closed for all 1,620 roles,
405 quartets, 135 inputs and 324 aliases. They retain 4,860 saved alignments,
19,440 candidate frames, 293,964,912 projected state observations and 3,453
unanchored candidate residues. All **388,800 historical rate cells** remain
unqualified. Twelve frames fail normalization; 4,848 pass the mean check
without gaining rate eligibility. The 216 unknown tip observations and
127 differences between separately logged conditional draws are explicit.
Those two logs do not establish a joint trajectory, and internal categories
remain unavailable in the historical output.

The completion archive binds 38,314 sources/artifacts and two original
producer/readback journals. These are output-accounting results, not mixing,
likelihood or biological acceptance.
[Full completed historical replay](../metadata/independent_short_sampler_replay_v3_completed_20261003.json).

Software checks retain twelve actual failed historical frames, two artificial
rate failures, two artificial unsuccessful roles, all 400 unknown-residue
support combinations, and full 1,620-role mocked producer/readback serialization.
Fourteen malformed-property cases and sixteen serialization cases are rejected.
The first new software invocation collided with a preexisting audit directory
before native execution. That failure is preserved; a separate fresh namespace
completed successfully with unchanged sources and no old output edits.
[Historical reader software proof](../metadata/independent_short_sampler_replay_v3_software_validation_20261003.json).

## Full corrected execution and resources

October 3, 15:20 UTC error recheck: all 1,620 corrected startup roles have
closed with zero unsuccessful roles. All 14,402 closure bindings were freshly
rehashed. The original V3 sampler is running under its declared limits; every
one of the 20 completed roles present at the audit snapshot passed read-only
native-output and serialized-array replay. This covers 60 joint frames and
3,883,978 ancestral residue/category pairs, including the unchanged strict
mean-one rate check. The installed original formatter still corrupts five of
twelve probe constants; native CJSON preserves all twelve. Full sampling and
posterior adequacy remain pending. The earlier startup/queued snapshots below
are retained as historical observations.
[Current native probe](../metadata/baliphy_logging_error_recheck_20261003_1520.json),
[saved-frame replay](../metadata/baliphy_current_joint_logger_error_recheck_20261003_1520.json),
[actual execution evidence](../metadata/baliphy_current_error_recheck_execution_20261003_1520.json).

The V5 startup grid is running native `--test` for every role under two CPUs,
32 GiB RAM, zero cgroup swap and two concurrent workers. Its producer, readback
and closure have new original handles. At 13:28 UTC, 216 startup checks had
passed with zero unsuccessful checkpoints; full closure remains pending.
[Original startup checkpoint](../metadata/baliphy_joint_logger_preflight_v5_execution_checkpoint_20261003_v3.json).

The V3 20-iteration joint sampler is queued behind four complete prerequisites:
this corrected startup, the closed historical native sampler, its closed
resource observer and the completed V3 historical replay. Its software gate
checks all 1,620 configurations, three retained corrected native fixtures,
full mocked workflow serialization with two failed roles, missing export
rejection and nine altered prerequisite cases. Historical-rate accounting is
mandatory; historical rates are never adopted by this new sampler.
[Joint workflow software gate](../metadata/baliphy_joint_sampler_software_validation_20261003_v3.json).

When admitted, the full grid uses sixteen CPUs, 200 GiB RAM, zero cgroup swap
and a 192-GiB reservation ledger. All 1,512 ordinary roles retain 12-GiB
address-space reservations; all 108 roles in the two risk families retain
48-GiB reservations. Reservations remain held through output audit. Every
failed/invalid role remains recorded, with no automatic native retry. Per-file
limits are 2 GiB; the free-disk guard is 256 GiB and planned output allowance
128 GiB. These are computational checks, not a justified longer posterior.
[Full sampler plan](../metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json)
and [queued original handles](../metadata/baliphy_joint_sampler_qualification_v3_launches_20261003.json).

A new read-only observer now follows the exact V3 sampler controller and every
role. Its producer uses one CPU/2 GiB/no swap and a one-second requested poll
cadence; readback/closure use two CPUs/4 GiB/no swap. Full configuration and
mock telemetry serialization passed, including two failed roles, two missing
live readings and eight altered/restart cases. The shared procfs decoder's
149 earlier actual native observations are reused through unchanged source
hashes; they are not current V3 readings. At 13:28 UTC the new observer and
queued sampler handles were freshly verified, with zero current native roles
and no counters reset. Its full observation/readback/closure remain pending.
[Observer plan](../metadata/baliphy_joint_sampler_resource_observation_v3_plan_20261003.json)
and [original runtime checkpoint](../metadata/baliphy_joint_sampler_resource_observation_v3_execution_checkpoint_20261003_v1.json).

The bounded software stages used two CPUs, 16 GiB and no cgroup swap. Recorded
child CPU seconds/peak RSS are 68.32/831,438,848 bytes for the native logger,
67.23/792,551,424 bytes for startup software, 12.26/233,594,880 bytes for the
historical reader software and 89.75/311,742,464 bytes for the joint workflow.
The resource-monitor software check used 55.73 child CPU seconds and
146,673,664 bytes child peak RSS. These child measurements are separate from
manager/cgroup peaks. The exact
original wait commands, child exits and invocation journals are retained in
[the transport manifest](../metadata/baliphy_cjson_workflow_software_original_transports_20261003.json).

GPU prediction remains paused. Neither old allocation failure is causally
resolved by these short runs. Posterior adequacy, full-grid corrected saved
frames, accepted phylogenetic/root/model controls and all eight biological
aims remain incomplete. No new infrastructure charges were incurred.
