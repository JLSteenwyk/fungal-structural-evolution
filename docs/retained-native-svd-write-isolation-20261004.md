# First factor change spans the native SVD call

The deeper original V6 diagnostic completed and caught a real input mutation
across the native LAPACK dispatch inside SciPy's SVD. The last unchanged event
precedes `_decomp_svd.py` line162 (`u, s, v, info = gesXd(...)`); the next
line165 event sees a changed `pmsf_profile_profile` factor. Its original probe
uses `profile_guide`, a different tree. This narrows the operation substantially;
the exact native write instruction and the original unprotected failure's cause
are not yet established. No source library, production tolerance or fit is changed.

The first four probes agree and preserve their inputs. Probe five raises the
explicit input-mutation diagnostic stop. All35unattempted identities stay
explicit. The eight changed values at row8513/columns115–122 and retained
snapshot are byte-identical to V2's captured array; its SHA256 remains
`d01b47ab2f25e8be93298456db41e20e03ccc9f360ad89292d89f4475a5347b1`.
Original wait82836 exits0; wrapper/invocation3a6c4b2756014618ba2b0d45d71741de,
whole terminal payload and manager start/end are verified. A successful
stop-and-capture diagnostic is not acceptance of changed inputs. Consumed
exports are freshly checked; full original2.38million-source verification is
inherited from the closed exact-cohort producer.

A separate GDB diagnostic now protects the interior pages of all five source
factors read-only. First/last partial pages remain writable to preserve adjacent
allocator metadata. It runs the unchanged original probe and records native
backtraces, fault address, instruction and exact process identities if a protected
fault occurs. Protection changes execution context; no original job is restarted.
The diagnostic controller preserves missing/failed target dispositions and GDB
exit codes, including post-exit inspection failures. It does not retry a native
fault or alter the mathematical guards.

At08:13UTC, original wait40959, wrapper1567650, GDB1567656 and inferior1567666
are live under invocation761385dcacf74166907922d21d90ef30. Actual2CPU32GiB/
noSwap limits and all five read-only protection ranges are verified in procfs.
Five probes have agreed with unchanged inputs at that snapshot. Source code,
GDB binary and source references are pinned; AS24GiB/7200CPU/14400wall/file256MiB
and2GiBoutput are upper diagnostic allowances, not a calibrated ETA.

Evidence:

- [Completed V6 capture](../metadata/retained_timing_ordered_probe_diagnostic_20261004_v6.json) and [original terminal proof](../metadata/retained_timing_ordered_probe_diagnostic_transport_20261004_v6.json).
- [First native SVD transition and snapshot binding](../metadata/retained_timing_probe_native_svd_transition_20261004_v6.json).
- [Protected debugger plan](../metadata/retained_factor_protection_debugger_plan_20261004_v1.json) and [exact live protection/identity check](../metadata/retained_factor_protection_debugger_checkpoint_20261004_v1.json).
- [Protected target](../scripts/diagnose_retained_timing_protected_factors_v7.py) and [failure-preserving debugger controller](../scripts/run_retained_factor_protection_debugger_v1.py).

Next: poll original40959, retain its native fault/completed-target evidence and
all probe dispositions, and verify the original terminal/source/artifact chain.
Validate any concrete backend repair independently before a fresh full retained
workflow. Weighted numerical work continues750/4340cohorts at08:11 with no
recorded failures. No real full weighted biological fit, adequate ancestral
posterior, GPU prediction resume or completed evolutionary aim.
