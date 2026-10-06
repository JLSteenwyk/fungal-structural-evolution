# Full domain extraction completion and independent atom readback

The refreshed coordinate extraction completed all **2,454,565 candidate
intervals from 976,357 original AFDB models**, in 977 shards. Every interval
was exported; the producer recorded zero rejected intervals and zero intervals
with missing backbone atoms. These are producer outcomes pending the separate
full atom reader, not accepted biological boundaries or homologous families.
All four annotation-policy alternatives, both alignment and envelope boundaries
and all 2,705,564 boundary associations remain unchanged.

The original four-worker extraction took **35,718.94 seconds (9 h 55 min)**.
Its archives contain **245,539,491,840 bytes**. The producer checked its complete
3,972-file binding map before returning. Coordinates and archives remain outside
Git; their paths and expected digests are retained in
`results/domains/full-domain-native-closure-20261006-v1/producer-source-bindings.json`.

## Original completion evidence

The original unit is
`fungal-completed-domain-coordinates-20261005-v1.service`, invocation
`5ded6cb02d854946ab9215ae72d2f9e4`. Its original wrapper PID 314168 and creation
time 1791223009.43 match the launch descriptor and the complete saved wrapper
initial and terminal journal messages. The execution receipt reports actual
native exit zero, no timeout, the complete bound producer receipt and unchanged
four-CPU/32-GiB/no-swap resource configuration. The original manager start and
CPU-accounting completion records are also present.

Original tool session **64749** is no longer available. A fresh poll returned
`Unknown process id 64749`; its API terminal exit code remains **null**. Native
completion is established by the saved execution and original journal, without
claiming recovery of that missing tool terminal. Extraction was not repeated.

The closure checker verifies the actual original records and rejects eleven
mutations of native exit, timeout, wrapper identity, invocation, child command,
API availability/zero, manager completion, wrapper terminal, shard scope and
global counts. The software check completed in its first tool return: no session
handle was assigned, and none is invented. Its exact return and native execution
are retained. The separate native-completion review, original session **92757**,
returned actual exit zero.

The first review-transport check, session **40520**, failed because the generic
checker required a CPU-accounting completion log that systemd omitted for the
818-ms review. This failure and its original records remain intact. The new
short-review transport checker verifies the actual original manager success,
zero exit and runtime from the original wait, alongside both complete wrapper
messages. It records the optional CPU-accounting count as zero. Its own original
session **7626** also returned actual exit zero. This correction changes review
bookkeeping; it does not modify the extraction or reconstruct session 64749.

Evidence:

- [Original producer result](../metadata/completed_domain_coordinates_20261005_v1.json)
  and [execution](../metadata/completed_domain_coordinates_execution_20261005_v1.json).
- [Native closure](../metadata/completed_domain_native_closure_20261006_v1.json)
  and [eleven rejection controls](../metadata/completed_domain_native_closure_software_20261006_v1.json).
- [Preserved transport failure](../metadata/completed_domain_native_closure_transport_failure_20261006_v1.json)
  and [successful review transport](../metadata/completed_domain_native_closure_transport_20261006_v2.json).

## Full independent reader complete; retained reporting review

The unchanged, previously qualified
`scripts/readback_domain_coordinate_archives.py` ran across the **entire
manifest**, behind an adapter requiring the closed original native execution
and the actual completed review transport. The original reader tool session was
**5441**; unit `fungal-full-domain-atom-readback-20261006-v1.service`, invocation
`07c83fa9e4b54a50800c1aebc9a6c2d8`, wrapper PID **749172**, creation time
**1791259825.33**. Its command, plan hash and actual kernel limits are recorded
in the [launch descriptor](../metadata/completed_domain_atom_readback_launch_20261006_v1.json).

The reader first reconstructs the full interval/job/catalog membership, with
no missing or duplicate intervals, model substitutions or foreign jobs. It then
freshly hashes every shard artifact and every original CIF used for export,
and checks all archive members, fragment sequences, atom/residue/element/chain
identities, coordinates, occupancy, confidence and confidence summaries.
The unchanged absolute tolerances are **0.000501** for coordinates, **0.005001**
for occupancy/confidence, and **1e-10** for confidence summaries. The separate
reader uses the same CIF lexical parser; recorded rejection reasons are retained
without claiming independent causal adjudication. The earlier fixture and its
coordinate-mutation rejection remain bound to the exact unchanged reader source.

Resources were estimated before launch: **eight CPU workers, 64 GiB memory,
zero swap, 56 GiB address-space cap per process, one BLAS thread**, a seven-day
wall cap and five-day per-process CPU cap. Reader metadata has a 1 GiB planning
allowance and 256 MiB per-file cap; minimum free space is 1 TiB. Launch observed
631.6 GiB available RAM and 9,758.7 GiB free disk. The **4–168 hour** planning
range is uncalibrated, not a finish ETA. The producer's runtime is context for
planning; reading and atom comparison have different costs. Source CIF I/O
volume is not yet estimated. No GPU, new prediction or paid resource is used.

The reader completed all **977 shards**, checking **976,357 models** and
**2,454,565 intervals**. Its retained raw receipt reports 2,454,565 exported
intervals, zero rejected intervals and zero intervals with missing backbone
atoms. The raw atom checks freshly hashed every source artifact and used CIF,
then checked every exported atom and its sequence, identity, coordinate,
occupancy, confidence and confidence summary against source arrays.

The original adapter then failed with its actual exit code one because its
strict dictionary equality check expected an explicit zero-valued `rejected`
key, while Python's `Counter` omitted it. This did not invalidate or rerun the
saved raw atom checks. The failed adapter execution, raw receipt and source
files remain unchanged.

A separate retained-output review completed against that exact terminal
execution and raw receipt. Its count handling passes every one of the 977
original producer shard records and eleven rejection controls. It accepts the
omitted `rejected` key only when the producer explicitly recorded zero, and
rejects every other missing, foreign, changed, noninteger or negative counter.
It repeats no atom checks. The review verifies the raw full-readback result;
it does not provide confidence/PAE calibration, biological boundary or homology
acceptance, or evolutionary inference.

[Observed reporting mismatch](../metadata/domain_atom_readback_reporting_defect_20261006_v1.json)
and [full count-scope software controls](../metadata/domain_atom_readback_count_scope_software_20261006_v2.json)
preserve this distinction. The [retained-output review](../metadata/completed_domain_atom_readback_count_review_20261006_v2.json)
records the accepted zero-counter scope.

The [reader plan](../metadata/completed_domain_atom_readback_plan_20261006_v1.json),
[resource estimate](../metadata/completed_domain_atom_readback_resources_20261006_v1.json)
and [original launch checkpoint](../metadata/completed_domain_atom_readback_checkpoint_20261006_goal_launch_v1.json)
and [partial atom-check checkpoint](../metadata/completed_domain_atom_readback_checkpoint_20261006_goal_progress_v2.json)
retain the complete scope and exact source digests.

```bash
python scripts/run_completed_domain_coordinate_atom_readback_v1.py \
  --plan metadata/completed_domain_atom_readback_plan_20261006_v1.json \
  --receipt metadata/completed_domain_atom_readback_20261006_v1.json
```

Reproduction requires the recorded original evidence and all pinned data,
the bounded wrapper/cgroup configuration, and fresh output/receipt paths.
Existing output paths are deliberately refused. Poll the original session or
observe its exact process/invocation; an observation failure is not permission
to restart it.

The full 501-fungal-entry plus 25-outgroup project remains active. Complete
prediction coverage, accepted phylogenetic/homology/copy/dating frameworks,
calibrated sequence–structure branch effects, adequate ancestral posterior
uncertainty and all eight evolutionary aims remain unfinished. GPU prediction
remains paused; this stage validates existing coordinates.
