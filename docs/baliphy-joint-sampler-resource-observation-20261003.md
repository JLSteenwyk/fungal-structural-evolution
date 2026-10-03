# Full joint-sampler resource observation

A separate read-only monitor now follows the captured original controller
for all **1,620 joint-logging qualification roles, 405 quartets, 135 inputs
and 324 aliases**. It starts while that controller is queued, so it can
observe native attempts when the four prerequisites close. It neither
launches nor signals inference and changes no native, GPU or cgroup limits.

The [immutable observation plan](../metadata/baliphy_joint_sampler_resource_observation_plan_20261003.json)
binds 1,437 sources, including the queued sampler plan, exact original launch
record, qualified observation code and verified software execution. At
10:26 UTC the observer/readback/closure handles 2978138/2978142/2978147 and
the exact original joint controller 2967091 were verified with their live
limits. The joint-stage controller had no native attempts or observations;
its small waiting-controller memory reading is not a native memory result.
[Original execution checkpoint](../metadata/baliphy_joint_sampler_resource_observation_execution_checkpoint_20261003_v1.json).

## Measurements and accounting

The unchanged, previously qualified procfs decoder checks every native
PID, creation time, argument list, configuration, command record, cgroup
and enforced address-space/CPU/file limit. It rechecks identity around
reads. Empty or transitioning argument lists, wrappers before native exec,
disappearing processes and missing fast-attempt readings remain unavailable
observations. A collected unit or missing process alone is not success;
the original controller's actual invocation-linked completion/failure
journals must be checked.

Per-native snapshots retain Linux-reported VmPeak, VmSize, VmHWM, VmRSS,
VmSwap and CPU time. Cgroup snapshots retain charged memory/current/peak,
swap and memory-limit/OOM event counters. Values are sampled sequentially,
not as a simultaneous concurrency census. Native RSS and charged group
memory have different scopes. No reset or configuration write occurs.

The monitor scans every role and records all attempt identities at original
controller termination, including attempts too fast to be observed live.
Failures, unstarted roles, missing receipts and missing observations remain
explicit. Peak fields remain null for unobserved roles; they are not
replaced by limits or inferred from successful exits. Polling gaps and
approximate high-water readings cannot establish exact final per-native
peaks or longer-chain memory requirements.

The serialized reader streams the complete event history, checks sequences,
limits, actual attempt metadata and original termination evidence, reconstructs
all 1,620 role records and compares every export/summary. Complete observer
source/artifact and two-original-journal closure are required. This is
observational closure, distinct from the sampler's own native output/readback
closure, likelihood/mixing qualification and biological acceptance.

## Software and resource qualification

The [software gate](../metadata/baliphy_joint_sampler_resource_observation_software_validation_20261003_v1.json)
checked every joint-stage configuration and reused the unchanged shared
procfs implementation's retained native fixture qualification. It passively
rechecked four actual currently running historical native workers. Those
are historical workers, not new joint-stage measurements; no inference was
launched. Original joint controller identity and its queued limits were
checked separately.

Full producer/reader serialization used all 1,620 real joint configurations
with artificial attempts, receipts and snapshots and explicitly mocked
capture/termination journals. Two artificial failures and two missing live
observations were retained. Eight restart/omission/alias/count/invented-peak/
scientific-claim cases were rejected. Nine shared decoder negative cases
remain bound through the prior qualified gate. These fixtures do not prove
full future native capture or actual observer closure.

Direct software execution completed with exit zero; all 52 bindings and
the captured log were verified. It used 54.74 seconds elapsed, 55.01 self
CPU seconds and 141,770,752 bytes self peak RSS. Its command configured
6 GiB address space, 150 CPU seconds, 180 seconds wall timeout, 256 MiB per
file and one-thread BLAS. Self memory is not charged group memory or a
native final peak.
[Execution evidence](../metadata/baliphy_joint_sampler_resource_software_execution_20261003.json).

The live observer has one CPU, 2 GiB RAM and no swap; readback and closure
each have two CPUs, 4 GiB and no swap. A one-second requested cadence
includes scheduling/processing delays, which are recorded. Planning allows
16 GiB output and requires 64 GiB free disk. The uncalibrated 0.1–4-day
planning sensitivity is not a finish ETA or enforced total-duration/output
limit. At preparation, 300.58 GiB RAM and 10,366.99 GiB disk were available.
No paid resources are provisioned.

Full joint sampler and observer native/output/provenance closure, longer
resource planning and adequate posterior ensembles remain open. Historical
allocation failures and model/root/likelihood/predictor/calibration concerns
remain unresolved. All eight biological aims remain incomplete. GPU
structure prediction stays paused.

Inspect original handles instead of relaunching:

```bash
python scripts/record_baliphy_joint_sampler_resource_observation_checkpoint.py \
  --output NEW_JOINT_RESOURCE_CHECKPOINT.json
```

Software checks use `scripts/check_baliphy_joint_sampler_resource_observation.py`
with a new audit root and receipt and the stated caps. Its passive current
worker check requires the original historical workers to be live; retained
past measurements are provenance evidence rather than new live readings.
Preparation uses
`scripts/prepare_and_launch_baliphy_joint_sampler_resource_observation.py`;
preserve existing plans, original handles and attempt/event roots.
