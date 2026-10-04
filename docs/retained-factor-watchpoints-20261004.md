# Preserve nonrecurrence and isolate an unprotected native write

The original retained timing failure remains unresolved. The unprotected V6
trace detected eight changed factor entries immediately after a native SVD
dispatch, retaining the changed array and all 35 unattempted identities.
This identifies a call spanning mutation, not the exact native write instruction
or an established library defect. See the
[original transition evidence](retained-native-svd-write-isolation-20261004.md).

A separate protected-memory run has now completed all 40 original selected
cohort-11 probes. All declared numerical points agree and all input hashes are
preserved. Original wait 40959 exited zero, with exact wrapper identity,
complete initial/terminal payloads, manager start/end and freshly checked
source/artifact hashes. The
[closed debugger receipt](../metadata/retained_factor_protection_debugger_20261004_v1.json)
and [terminal proof](../metadata/retained_factor_protection_debugger_transport_20261004_v1.json)
preserve this outcome. No SIGSEGV or protected-page write was observed.

Read-only protection applies to complete interior pages of the five factor
buffers; edge pages stay writable. This changes execution/allocation context.
Nonrecurrence under protection is not repair and does not invalidate the
earlier unprotected corruption. GDB's post-normal-exit inspection commands
print errors because the inferior has exited; those exact logs are retained.
There is no native fault address, fault backtrace or implicated instruction
from this successful protected run.

A new unprotected hardware-watchpoint diagnostic now runs the preserved
original probe, without read-only page protection or numerical changes.
Two aligned eight-byte words, row 8513/column 116 in `pmsf_profile_profile`
and `pmsf_mafft_profile`, cover the interior of the previously changed spans.
The inferior writes its exact live buffer addresses, raises SIGSTOP, then GDB
arms two hardware watchpoints and continues. The first observed write stops
execution and retains old/new values, native stack, registers and disassembly.
The watchpoints do not cover every byte of all five buffers. Debugger/SIGSTOP
still change execution context, so nonrecurrence would remain inconclusive.

Before launching the actual factor diagnostic, an
[artificial positive control](../metadata/retained_factor_watchpoint_control_20261004_v1.json)
successfully installed two hardware word watchpoints and caught a deliberate
native `memset` write on private synthetic words. Its
[original terminal proof](../metadata/retained_factor_watchpoint_control_transport_20261004_v1.json)
verifies actual wait 58131, exit zero, whole wrapper payloads and manager
start/end. This qualifies capture mechanics; it is not observed source-factor
corruption or a biological result.

The [new actual plan](../metadata/retained_factor_watchpoint_debugger_plan_20261004_v1.json)
binds 78 sources and 49 project modules, unchanged cohort exports and both
earlier diagnostic/reference closures. Two CPUs, 32 GiB/no swap, 24 GiB
address space, 7,200 CPU seconds, 14,400 wall seconds and 256 MiB per-file
limits bound this separate diagnostic; output planning is 2 GiB. No original
failed stage is restarted. Native libraries, priors, optimizer settings and
production numerical tolerances remain unchanged.

At 08:54 UTC, the
[exact original live checkpoint](../metadata/retained_factor_watchpoint_debugger_checkpoint_20261004_v1.json)
verifies wait 91413's live wrapper 1607519, GDB 1607528 and inferior 1607538,
actual caps, one BLAS thread and both watched words mapped writable. The
first two probes agree with all input hashes preserved. No final receipt or
native write instruction is established at this checkpoint.

Scripts:

```text
scripts/diagnose_retained_timing_watched_factors_v8.py
scripts/run_retained_factor_watchpoint_debugger_v1.py
scripts/retained_factor_watchpoint_positive_control_v1.py
scripts/record_retained_factor_watchpoint_checkpoint_v1.py
scripts/record_capped_watchpoint_debugger_transport_v1.py
```

The protected receipt transport uses a new V2 helper to read the correct
debugger output root; V1 was executed only for the preparer branch and remains
preserved. No prior executed source was edited.

Original guard failure, partial timing outputs and all failed downstream
controllers remain immutable. No full-data weighted biological fit has run;
weighted numerical work continues independently. Finding an operation or write
instruction would still require a separately verified repair and complete
source/numerical/timing gates before biological fitting.
