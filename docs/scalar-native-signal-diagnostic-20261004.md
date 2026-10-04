# Preserved V6 native sampler failures — October 4

The 11:27 original sampler checkpoint has 1,174 of 1,620 recorded outcomes:
1,152 finite output-integrity checks, 10 explicit special-value reviews and
12 unsuccessful attempts. Full sampler/telemetry readback and closure remain
pending. Twenty iterations do not establish posterior uncertainty.

All twelve presently saved failures belong to **OG0000972**, one effective
input group shared by two configuration aliases. Every chain role, across
broad, centered and package priors, exits **-11 (SIGSEGV)**. Native receipts,
original PID/create/command identities and saved artifacts are independently
checked. Stdout reaches “Beginning MCMC computations”; stderr is empty and
no allocation warning or bad_alloc message was recorded. Failed roles admit
no ancestral arrays and remain in full accounting.

A bounded original telemetry prefix gives 1,185 matching observations for
these exact twelve identities. Maximum sampled VmPeak is 35.61 GiB and VmHWM
35.23 GiB under their original 48 GiB address-space limits. Sampled values can
miss later peaks and do not prove sufficient memory or a cause. The extracted
snapshots remain separate immutable evidence; the growing original telemetry
file is not described as having a stable full-file hash.

- [Original sampler checkpoint](../metadata/baliphy_scalar_v6_sampler_execution_checkpoint_20261004_goal_1128.json)
- [Full saved failure/identity/artifact recheck](../metadata/baliphy_scalar_v6_failed_attempts_recheck_20261004_v1.json)
- [Matched original telemetry extract](../metadata/baliphy_scalar_v6_failed_attempts_telemetry_extract_20261004_v1.jsonl)

Two private actual GDB controls qualify signal-stack capture and normal exit:
an artificial ctypes null-pointer access stops at SIGSEGV with a native stack;
a normal target exits without invalid register inspection. Original wait
77247, whole wrapper messages, source hashes and manager journals are closed.
[Software transport](../metadata/scalar_native_signal_debugger_software_transport_20261004_v1.json).

A separate debugger inferior uses the first failed broad-prior role's original
program, seed, iteration count and native limits: **48 GiB AS, 5,674 CPU seconds,
2 GiB per file**. It writes under a fresh diagnostic output root. GDB and its
wrapper share two CPUs, 64 GiB memory, no swap and one BLAS thread. Original
failed attempts and all remaining full-grid roles are unchanged. The observed
failed-role median was 109 seconds; the prospective 2–15 minute debugger range
is uncalibrated, with a two-hour native observation limit. No core dump is
requested, no infrastructure is provisioned and no charges are incurred.

The driver captures the signal, PC, siginfo, backtrace, registers, instructions,
shared libraries and `/proc` status/limits/maps if the native target stops.
Debugger execution changes context; neither failure nor nonrecurrence alone
establishes a repair. A captured diagnostic stop is not a completed sampler.

- [Frozen diagnostic plan](../metadata/scalar_native_signal_debugger_plan_20261004_v1.json)
- [Exact original wrapper and observed native limits](../metadata/scalar_native_signal_debugger_checkpoint_20261004_v1.json)
- Original tool wait: **19850**; invocation `f1261b24602246d899fd02bc16b84c9d`.
- Output root: `results/ancestral/scalar-native-signal-debugger-20261004-v1/`.

The original wait 19850 finishes with exit zero at 11:39 after capturing a
SIGSEGV, not a completed native sampler. Original wrapper/manager transport and
all diagnostic artifacts are verified. The PC is
`reg_heap::incremental_evaluate1_changeable_(int)+56`, at a call instruction.
The fault address equals `rsp-8`, 216 bytes below an 8 MiB stack mapping. Native
VmStk and its soft limit both equal 8 MiB; the hard limit is unlimited. Repeated
evaluator frames support stack exhaustion in this separate diagnostic. This
does not prove that all twelve original crashes share the cause.
[Verified instruction/mapping/limits evidence](../metadata/scalar_native_stack_fault_20261004_v1.json).
[Original terminal and journal proof](../metadata/scalar_native_signal_debugger_transport_20261004_v1.json).

A separate V2 driver also records actual native exit codes. Real private
signal/normal-exit controls pass with a 64 MiB process-local stack, original
wait 92249 and full journal proof. A controlled diagnostic now uses that stack
soft limit, retaining unlimited hard stack and all original AS/CPU/file caps,
program/input/prior/seed. Production limits and all original output bytes
remain unchanged. Native command identity and actual process limits are checked
directly. Original wait **99691**, invocation
`eadee75fd20e4b51ab6d9ab87c2a7af1`, remains pending; no validated fix or posterior
acceptance is inferred.

- [V2 actual controls](../metadata/scalar_native_signal_debugger_software_transport_20261004_v2.json)
- [Frozen controlled comparison](../metadata/scalar_native_stack_comparison_plan_20261004_v2.json)
- [Original comparison identities and limits](../metadata/scalar_native_stack_comparison_checkpoint_20261004_v2.json)

Poll the same original wait for terminal evidence; observation timeout is not
grounds for restarting. No source, prior, tolerance, seed or process cap is
changed in the production sampler. Full independent output integrity and
adequate posterior sampling remain required, alongside all eight aims.


## Controlled 64 MiB stack outcome, 13:28 UTC

The unchanged-input/prior/seed/20-iteration comparison ends with native
SIGKILL. Original debugger wait 99691 exits zero, but its native child never
completes the planned horizon. The original callback's null code/signal and
normal-exit inspection label do not prove normal native termination. Original
GDB text is explicit; every artifact and whole wrapper/manager record is
verified. Ten saved scalar rows, iterations 0 through 9, independently compare
all 430 mapped values. No posterior arrays are adopted from this partial run.

The 5,674 CPU-second hard cap and process-local 64 MiB soft/unlimited-hard stack
remain unchanged. Prior live CPU consumption approached the cap before the
child disappears; cap attribution is consistent with observations rather
than separately captured kernel evidence. The stack increase alone is not a
validated correction for the original 24 SIGSEGV failures. No automatic retry,
production policy or machine-wide limit is changed.
[Full terminal review](../metadata/scalar_native_stack_comparison_terminal_review_20261004_v3.json).
