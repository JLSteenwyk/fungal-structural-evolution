# Full future joint-logger startup validation

Native startup validation is running for **all 1,620 future roles, 405
quartets, 135 effective inputs and 324 original aliases**. Every command uses
`--test --log-format json`; no MCMC iterations are requested. This qualifies
the prepared source programs and initial model representation before a
future sampling horizon. It does not test saved full-grid joint logging
frames, adequate mixing or ancestral structural predictions.

## Inputs, invariants and accounting

The [prepared future sources](../metadata/baliphy_joint_node_logger_future_models_20261003_v4.json)
are reused without editing. All future seeds are unique and disjoint from
the 3,243 earlier real/software seeds. Each startup job preserves its old
short-sampler source identity, aliases, alignment/tree hashes, family, prior
and role. Its new program reverses exactly to the old reference-initialized
program through the three qualified logger edits. Startup seeds are not
continued posterior chains or additional independent biological replicates.

For each attempt, the producer checks configuration/command/process identity
and all artifacts. It verifies the source inverse, rooted clade and branch
correspondence, observed reference alignment/homology, fixed-tip versus
free-ancestor representation and degree-aware indel density. It accepts one
initial-model JSON object with only the exact recognized native timing footer.
All failed/invalid native startups remain explicit. Exit zero alone is not
accepted model initialization.

Readback repeats every native and serialized check for every role, verifies
all exports and hashes, and compares the complete summary. The reader shares
the decoder; it is not a third independent likelihood implementation. Complete
source/artifact hashes and both original producer/reader completion journals
are required for accounting closure. Completed stages refuse restart; existing
attempts are never silently replaced or selected among multiple attempts.

## Software qualification

The [full startup software qualification](../metadata/baliphy_joint_logger_preflight_software_validation_20261003_v1.json)
passed complete real-role/source/configuration generation for all 1,620 roles
and 405 programs. Eleven altered designs were rejected: missing/duplicate
roles or source identities, duplicate/reused/invalid seeds, changed roles,
aliases, priors or source hashes, and invented execution claims.

Three capped native synthetic `--test` fixtures, one per prior, validated
reference homology, density, rooted tree and fixed-tip/free-ancestor behavior
using the actual corrected programs. Full producer/reader serialization then
used every real role's metadata with explicitly mocked native execution and
source admission. Two artificial failures and two unresolved quartets remain
retained; completed restarts, missing checkpoints, changed seeds, scientific
acceptance and false posterior claims were rejected. These mocks do not prove
native startup success for the full fungal case grid. The fixtures are
computational QC, not a separate fungal pilot.

The original software controller 2931768 and its actual invocation-linked
terminal/resource journal were verified, together with 41 source bindings
and its captured two-CPU/16-GiB/no-swap group limits. It used 73.989164 CPU
seconds, 73,945,088 bytes peak charged group memory and zero swap. Each of its
three native fixtures had 8 GiB address space, 120 CPU seconds, 180 wall
seconds and a 64 MiB file cap. Planning allowed 2 GiB output and required
64 GiB free disk. No posterior sampling, GPU or charges occurred.
[Verified software execution](../metadata/baliphy_joint_logger_preflight_software_execution_20261003.json).

## Full production startup and resource plan

The [immutable plan](../metadata/baliphy_joint_logger_preflight_plan_20261003.json)
binds 1,432 sources, including the actual software execution proof. Its jobs
are under `results/ancestral/joint-logger-startup-inputs-20261003-v1` and its
new output root is `results/ancestral/full-baliphy-joint-logger-preflight-20261003-v1`.
The preparation record and all older programs, sources, jobs and attempts
remain unchanged. The historical short sampler and both output readers keep
their original schedules and resource caps.

The [launch inventory](../metadata/baliphy_joint_logger_preflight_launches_20261003.json)
records controllers 2935520/2935524/2935529 for producer/readback/closure.
Each has two CPUs, 32 GiB RAM, no swap and one-thread BLAS settings. The
producer runs two native workers, each capped at 12 GiB address space,
600 CPU seconds, 900 wall seconds and 256 MiB per file. Eight GiB of group
headroom remains beyond the two address-space caps. Planning allows 16 GiB
output and requires 64 GiB free disk; output estimates and enforced limits
are distinct. No paid resources are provisioned.

The complete earlier startup census measured 12,554.4245 elapsed worker
seconds across 1,620 roles, about 3.49 worker hours. Planning sensitivity is
6,277–50,218 worker seconds for the new sources/seeds. This excludes a
calibrated assessment of the current shared host load and is not a finish
ETA or memory guarantee. All configured wall timeouts total 405 worker
hours as an upper allowance, not expected runtime. At preparation, 442.24
GiB RAM and 10,370.62 GiB disk were available.

At 09:18 UTC all three original handles, live caps and 1,432 pins were
verified; 12 startup roles had valid unclosed dispositions across all priors.
Full production startup/readback/closure remains pending.
[Original execution checkpoint](../metadata/baliphy_joint_logger_preflight_execution_checkpoint_20261003_v1.json).
At 09:20 UTC the separate original short sampler had 1,068 successful
unclosed role checks, while its observation-aware replay remained queued
with zero role checks.
[Short-sampler/replay handles](../metadata/independent_short_sampler_replay_v2_execution_checkpoint_20261003_v4.json).
Its [resource observer](../metadata/baliphy_sampler_resource_observation_execution_checkpoint_20261003_v7.json)
was reverified without resets or configuration changes.

The [joint logger software checks](baliphy-joint-node-logger-20261003.md)
demonstrate same-record sequence/state consistency in their capped native
fixtures. This startup-only gate cannot extend that result to full native
sampling output. Full-grid joint-frame qualification, sampler/resource and
historical reader closures, allocation failures, longer-chain resources and
mixing, likelihood/root/model/predictor controls, calibration, phylogeny,
reconciliation, dating and all eight biological aims remain open. GPU
structure prediction stays paused; no longer posterior horizon is launched.

## Reproduction

```bash
python scripts/record_baliphy_joint_logger_preflight_checkpoint.py \
  --output NEW_JOINT_STARTUP_CHECKPOINT.json
```

Inspect original handles instead of relaunching frozen stages. Software
checks can use new audit/receipt paths with
`scripts/check_baliphy_joint_logger_preflight.py`; preserve failed attempts
and the documented group/native caps. Source preparation uses
`scripts/prepare_baliphy_joint_logger_preflight.py` and requires the passed
full software gate. Original software terminal proof and resource records
were added to the unlaunched plan before its first launch and hash binding.
The [existing Python environment](../environments/independent-short-sampler-replay-20261003.yml)
records dependencies; native executable/API/source hashes are bound in the plan.
