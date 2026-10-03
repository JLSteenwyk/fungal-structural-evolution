# Native ancestral sampler stack exhaustion and scoped correction

The complete corrected joint-logger sampler has closed all **1,620 roles,
405 quartets, 135 inputs and 324 configuration aliases**. Independent output
readback accepted 1,596 short attempts and retained 24 native SIGSEGV failures:
399 complete and six unresolved quartets. All failures belong to two
622-protein whole-alignment inputs in `OG0000972`, with branch floors 1e-7
and 1e-9 at resolution zero. Both retain all three priors and four chain roles.

The [sampler closure](../metadata/baliphy_joint_sampler_qualification_v3_completed_20261003.json)
binds 84,071 source/artifact hashes and both original producer/readback
journals. Its 4,788 joint frames contain 54,527,895 ancestral residue/category
pairs. These 20-iteration checks are **not adequate ancestral posteriors**.
The installed numeric formatter remains unchanged; the qualified V5 CJSON
logger bypasses it. Historical malformed rate records and two earlier
allocation failures remain preserved.

## Causal diagnostics

One original broad-prior chain from each failed input was examined in new
debugger namespaces. Commands kept the original model, alignment, tree, seed
and 48-GiB address-space/2-GiB file/native CPU caps. Diagnostics used
`--iterations 0` to isolate initialization and the first saved record;
core files were disabled. Diagnostic samples never enter posterior ensembles.

Both default-stack diagnostics reproduced SIGSEGV in the native evaluator.
The mapped stack was 8 MiB; the stack pointer was 144 and 176 bytes below
its lower boundary. Adding only `--stack=67108864` to each native resource
command allowed both processes to exit normally and write the initial
alignment/joint record. Independent decoding checked all 622 observed tips,
native ancestors, strict mean-one rates, states, residue coordinates and four
candidate clades. This establishes the initialization mechanism in two
representatives, not full sampling success in every role.
[Two-input diagnosis and complete failure accounting](../metadata/baliphy_stack_followup_diagnosis_20261003_v1.json).

The first debugger preparation could not set a nested file limit and never
executed BAli-Phy; it remains preserved. A fresh namespace corrected that
wrapper limit. Debugger exit zero is distinguished from its inferior's
SIGSEGV using explicit debugger events and retained native artifacts. The
second-input diagnostics also retain original controller identities,
invocation-linked journals and observed terminal tool returns.

The [complete original resource observer](../metadata/baliphy_joint_sampler_resource_observation_v3_completed_20261003.json)
closed all 1,620 attempt identities, with 88,437 live observations and 545
unavailable observations retained. It binds 7,933 hashes and both original
observer journals. Maximum sampled cgroup reported peak was 129,469,779,968
bytes (approximately 120.6 GiB). The full journal passed unchanged-limit,
zero-swap and no-observed-OOM/memory-limit-event checks. These samples are
not precise final native peaks or longer-chain memory guarantees.

## Fresh follow-up retaining the complete design

The [immutable plan](../metadata/baliphy_stack_followup_plan_20261003_v1.json)
contains exactly one fresh 20-iteration attempt for every failed role.
Programs, priors, reference initialization, alignments, trees, native
address-space/file/CPU/wall caps and candidate mappings are unchanged. Only
the scoped 64-MiB stack flag and a fresh seed differ. Seeds occupy a separate
namespace, checked against all original/current/historical ensemble seeds
and listed software fixture seeds. New identifiers/output roots prevent
continuation, concatenation, restart or replacement of original samples.

A full 1,620-role ledger references every original checkpoint by hash.
Successful originals stay immutable references. Each original failure keeps
its native receipt and disposition alongside the separate fresh outcome.
New failures remain unresolved without automatic retry. Complete quartet
accounting uses each applicable follow-up outcome while retaining the entire
original ledger. No taxa, priors, families or original settings are excluded.

Four workers hold 48-GiB FIFO reservations through execution and output audit,
sharing 192 GiB under four CPUs/200 GiB/no swap and single-thread BLAS.
Identity-bound native stack/memory and cgroup events are sampled throughout;
missing/transitioning readings stay explicit. Independent readback checks
native identity/configuration/artifacts, scalar rows, joint frames at
iterations 0/10/20, candidate arrays, reservations, resources and every
original/fresh role record. Missing exports cannot be recreated. Complete
source/artifact and both original new completion journals precede closure.

Readback and closure each have two CPUs/32 GiB/no swap. Preparation required
200 GiB available RAM and 256 GiB free disk; 311.9 GiB RAM and 10,338 GiB disk
were available. Output planning allows 32 GiB excluding referenced originals;
this is not a hard quota. Configured timeouts total 41.83 worker hours, an
upper allowance rather than an ETA. Initial logging does not calibrate full
sampler runtime. No new infrastructure or charges are involved.

## Qualification and current execution

The [software gate](../metadata/baliphy_stack_followup_software_validation_20261003_v1.json)
checks all original metadata and 24 generated follow-ups. Full producer/reader
serialization used explicit native, source-closure, cgroup and procfs mocks.
Three artificial fresh failures stay accounted, leaving 1,617 selected
successful and three unresolved software roles. Nine altered designs and ten
serialization/resource/restart cases were rejected, including coverage,
alias, prior, failure-class and accounting alterations; missing exports;
wrong stack limits; changed commands/identities; and OOM events. Copied
original arrays are software fixtures, not fresh fungal ancestral samples.

The actual software tool wait exited zero; exact original process and
invocation-linked completion evidence appear in the
[execution record](../metadata/baliphy_stack_followup_software_execution_20261003_v1.json).
Configured limits were two CPUs/8 GiB/no swap, 6 GiB address space, 480 CPU
seconds, 600 service wall seconds and 1 GiB per file. Launch settings are not
retrospective live procfs measurements; the manager's small completion-memory
value is not interpreted as native peak RSS.

Full original sampler/resource archives were freshly checked before launch.
Producer/readback/closure controllers are 3330319/3330323/3330328. At
**18:28 UTC**, all five original stage/dependency handles and caps were
verified. Four native workers had the expected 64-MiB stack; no fresh role
had completed yet.
[Original runtime checkpoint](../metadata/baliphy_stack_followup_execution_checkpoint_20261003_v1.json).

Recheck without relaunching:

```bash
python scripts/record_baliphy_stack_followup_checkpoint.py \
  --output NEW_STACK_FOLLOWUP_CHECKPOINT.json
```

Reproduce software QC in new paths with `scripts/check_baliphy_stack_followup.py`.
Native execution is `scripts/run_baliphy_stack_followup.py`; preparation and
launch are `scripts/prepare_and_launch_baliphy_stack_followup.py`. Qualified
sources, plans, attempts and counters stay immutable.

Full follow-up sampling/readback/closure, longer-run resources, adequate
posteriors, model/root/predictor qualification, accepted phylogenies,
reconciliation/dating and calibrated evolutionary effects remain pending.
All eight biological aims remain incomplete. GPU prediction stays paused.
This correction changes no global stack, SLURM, GPU or power default.
