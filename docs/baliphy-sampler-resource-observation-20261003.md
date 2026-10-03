# Resource observation of the full ancestral sampler

A separate read-only observer is recording the original full 1,620-role
short sampler. It tracks actual native process identities, enforced limits,
Linux memory readings, CPU use and whole-group memory events to inform the
resource plan for longer sampling. It launches no additional native inference
and changes no inference limits. Its complete serialized and provenance
closure remains pending while the original sampler runs.

## Scope and resource budget

The [immutable observer plan](../metadata/baliphy_sampler_resource_observation_plan_20261003.json)
binds 1,000 source/input/configuration hashes, including the original full
sampler plan and launch. It retains all 135 inputs, three priors, four roles,
405 quartets and 324 original aliases. Fast, failed, unstarted and unobserved
attempts remain explicit; none are replaced with invented memory values.

| Observer setting | Value |
| --- | --- |
| Producer CPU / memory / swap | 1 CPU / 2 GiB / zero |
| Reader and closure, each | 2 CPUs / 4 GiB / zero swap |
| Nominal poll interval | 1 second, plus processing and scheduling delay |
| Planned output allowance / minimum free disk | 16 GiB / 64 GiB |
| BLAS threads / new native jobs / GPU / new charges | 1 / 0 / none / none |
| Runtime | Depends on the original sampler; no calibrated finish ETA |

The planning scenario spans 0.1–4 days and is not a runtime guarantee. Actual
poll gaps are retained. The allowance is a storage estimate; the free-disk
guard is enforced. Streaming readback avoids retaining the entire journal in
RAM. Existing native roles retain their full 16-CPU/200-GiB/no-swap group,
192-GiB shared reservations and family-wide 12/48-GiB address-space limits.
The observer has no ability in its code path to reset counters, signal or
restart inference, change cgroup limits or alter GPU settings.

## Measurements and identity checks

Every discovered attempt binds its configuration, exact command, PID,
creation time, process-group identity and source hashes. A memory sample is
accepted only while the native identity, cgroup membership and enforced
address-space, CPU and file-size limits match, with identity checked before
and after the operating-system reads. The pre-exec wrapper is distinguished
from the native executable. Reused PIDs, empty/transitioning commands,
disappeared processes and zombies yield unavailable observations; a stable
unrelated command is rejected.

Per-process records retain VmPeak, VmSize, VmHWM, VmRSS, VmSwap, CPU user/system
time, and hashes of the raw status/limit readings. Linux RSS accounting can be
approximate. The last live sample can precede an attempt's final peak;
reported high-water readings are not exact final per-attempt memory bounds.
[Linux proc documentation](https://docs.kernel.org/filesystems/proc.html).

Whole-group readings retain memory.current, memory.peak, swap and memory.events.
The group includes descendants and kernel/cache charges, so its memory is
not a sum of native RSS. Peak readings are historical group counters; the
observer never resets them. [Linux cgroup v2 documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html).

Native readings are sequential, not simultaneous concurrency measurements.
The sampler's independent reservation ledger enforces admission and remains
the evidence for the 192-GiB reservation budget. Resource snapshots cannot
establish long-chain safety, convergence, a repaired allocator or an accepted
ancestral model.

## Complete accounting and software evidence

The producer records each polling event in observations.jsonl, preserving
actual gaps and unavailable observations. At original-controller termination,
it requires the exact original invocation-linked terminal journal and captures
every attempt descriptor, including attempts missed between polls. It writes
one disposition for every requested role. Full streaming replay binds native
configuration/command/process/receipt hashes, reconstructs live/missed counts
and reported maxima, and preserves failures. A successful original sampler
also requires all 1,620 attempt identities and receipts; a failed controller
retains explicit missing/unstarted roles.

The reader reconstructs every row, summary and source binding from the journal;
final observer closure requires both original observer completion journals and
the full source/artifact archive. Readback verifies saved-record integrity,
not the historical operating-system state again. Observational closure is
separate from the sampler's output closure and scientific qualification.

The [combined software evidence](../metadata/baliphy_sampler_resource_observation_combined_validation_20261003.json)
checks all 1,620 real role configurations, missing-observation accounting,
nine malformed/altered observation rejections and eight serialization/restart
rejections. One actual capped synthetic five-tip native fixture produced
149 verified live readings and exited successfully. Full producer/reader
serialization used explicitly artificial process/receipt/OS records, retaining
two failures and two roles without live observations. Those artificial records
are not production measurements or a biological pilot.

The first software observer failed its native-command identity assertion even
though the native fixture exited successfully. Its unexpected command was
not retained, so the precise triggering state is unproven. Its
[source and journal remain preserved](../metadata/baliphy_sampler_resource_observation_development_failure_20261003.json).
Before production launch, explicit empty/exec/exit handling was added and
tested alongside rejection of stable wrong commands. Both revised software
invocations have [actual original completion/resource records](../metadata/ancestral_qualification_completion_execution_checkpoint_20261003.json).
No production job was restarted by this repair.

## Observed execution and reproduction

At 07:13 UTC (03:13 EDT), a [fresh exact-handle observation](../metadata/baliphy_sampler_resource_observation_execution_checkpoint_20261003_0714.json)
verified all 1,000 pins, all three original observer handles and the original
sampler controller. The latest complete snapshot was sequence 1,313 and the
journal held 7,107,039 bytes. Original group memory.current was 6,174,265,344
bytes (about 5.75 GiB), memory.peak 11,863,592,960 bytes (about 11.05 GiB),
swap zero, with no recorded pressure/limit/OOM events. These are intermediate
group readings, not final native peaks or evidence that later roles will fit.

At 07:15 UTC, [19 sampler roles had successful unclosed output checks](../metadata/baliphy_reference_sampler_execution_checkpoint_20261003_v3.json).
Four exact native workers were live with 48-GiB address-space caps. Their
reservations use the entire 192-GiB budget; the maximum 16-worker setting does
not override the memory admission rule. The historical two failed chains
remain unresolved. Full sampler and observer closure, independent joint
site-property replay, longer disjoint-seed sampling and posterior/model/root
qualification remain open. GPU inference stays paused and all eight biological
aims remain incomplete.

```bash
python scripts/record_baliphy_sampler_resource_observation_checkpoint.py \
  --output NEW_RESOURCE_OBSERVATION_CHECKPOINT.json
python scripts/record_baliphy_reference_sampler_checkpoint.py \
  --output NEW_SAMPLER_CHECKPOINT.json
```

Use new output names. Do not rerun the producer/launcher while the original
observer is active or after its completed output exists. Launched source,
plans, native invocations and outputs remain immutable. The
[CPU environment](../environments/baliphy-sampler-resource-observation-20261003.yml)
records dependencies; native binary/API/model hashes remain bound by the
sampler plan. Large native artifacts and observation journals remain outside
Git with provenance in plans and completed archives.
