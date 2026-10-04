# Traced native SVD watchpoints — October 4

The earlier retained-input covariance error remains unresolved. The V8 target
completed 40 unchanged probes with no factor mutation or watchpoint hit, but
its original GDB/wrapper wait exited 1 because register inspection ran after
the inferior had exited. That original failure is preserved; it is not a
successful debugger completion or a repair.

A separate V2 debugger driver checks whether the inferior still exists before
requesting registers. Two actual GDB controls pass: a deliberate native write
to private words is caught, and normal exit avoids post-exit inspection. Both
controls and their original terminal/journal transport are verified.
[Control evidence](../metadata/retained_factor_watchpoint_debugger_v2_software_transport_20261004_v1.json).

The fresh V9 diagnostic combines two hardware watchpoints on the previously
affected source-factor words with Python line tracing around the original
probe. Source words remain mapped writable; no mprotect or model/guard changes
are introduced. Actual pre-call metadata records operand shape, dtype, strides,
pointer, native driver/callable, workspace and source aliasing. It does not
capture LAPACK's internal copied operands or work buffers.

Three actual small SciPy gesvd calls qualify metadata for C order, Fortran
order and direct source aliasing. Python emits repeated line events around
one call; repeated metadata events are not additional native calls. The first
checker wrongly required one line event and failed. Its original failure is
preserved; the separate V2 checker verifies identical repeated headers without
changing the diagnostic or numerical code.
[Metadata qualification](../metadata/retained_factor_v9_dispatch_software_transport_20261004_v2.json).

The actual original wait **35504**, invocation
`f9601cc5bc8b448683104fe9254f4c9f`, completed with exit **0** at 11:12 UTC.
All **40 of 40** ordered groups agreed and preserved their inputs; no native
source write was caught. All 46 target artifacts are freshly hashed, and
the whole original wrapper initial/terminal payloads and manager start/end
records match. The original exec initial tool payload was not retained; the
original terminal payload and wrapper initial message are preserved.
An observed dispatch used a float64 22,881-by-10 C-order operand,
gesvd workspace 22,911, and no alias with any of the five cached source factors.
Other calls have different dimensions; the 10:15 snapshot's final dispatch was
a 5-by-5 singular-value-only call.

- [Frozen plan](../metadata/retained_factor_watched_traced_debugger_plan_20261004_v1.json)
- [Exact original live identities, limits and trace observation](../metadata/retained_factor_watched_traced_debugger_checkpoint_20261004_goal_1015.json)
- [Completed original diagnostic and artifact transport](../metadata/retained_factor_watched_traced_debugger_transport_20261004_v1.json)
- [Fresh formatter and corrected-value check](../metadata/current_analysis_error_recheck_20261004_user_1122.json)
- Target outputs: `results/phylogeny/retained-timing-watched-traced-factor-probes-20261004-v9/`.
- Debugger outputs: `results/phylogeny/retained-factor-watched-traced-debugger-20261004-v1/`.

The diagnostic uses two CPUs, 32 GiB memory, zero swap and one BLAS thread,
with a two-hour native CPU and four-hour wall cap. Original terminal evidence
was obtained by polling the same wait; no diagnostic was restarted.
Forty nonrecurrences under this execution context do not prove that
the original unprotected corruption is fixed. No exact write instruction,
upstream library bug, ABI mismatch or installed-software repair is established.
