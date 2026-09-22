# Reconciliation stall diagnosis and traced recovery

The first full expanded reconciliation run is terminal and incomplete. The
native process returned zero after reporting a 120-second progress stall at
70,306 of 70,307 queued families. The execution controller rejected the run
because mandatory final outputs were absent; the MAFFT-guide stage did not
start. No final reconciliation result is accepted from this execution.

Inspection found 70,306 nonempty intermediate resolved-tree files. The sole
missing queued family is OG0000017, containing 55,981 proteins. These file
counts do not independently validate the intermediate trees. All 96,016
copied source input files still match their frozen manifest. Evidence is in
`metadata/expanded_reconciliation_v1_failure_review.json`; failed outputs
remain at `results/orthology/expanded-reconciliation-execution-v1/`.

The installed OrthoFinder 3.1.5 function `RunOrthologsParallel` defaults to
`STALL_TIMEOUT=120.0` and measures time since a task-completion message. A
single large unfinished task can therefore trigger it even if a worker is
still computing. The failed run did not record the worker's state at that
moment, so slow computation versus another stall is not yet established.
OG0000017 also retains its separate annotation/homology concern; longer runtime
does not resolve that scientific issue and the family is not silently omitted.

`scripts/run_orthofinder_with_progress.py` adds a run-specific wrapper around
the native launcher. It passes a 21,600-second no-task-completion allowance
and writes per-worker JSONL start/return/error events with family ID, process
ID, elapsed time and CPU time. Original native inference functions still
receive the same data and return their original results. The installed
package is unchanged. A recorded function return is not proof of valid
inference, especially because native code catches some internal exceptions;
mandatory outputs and independent output readback remain required.

The actual four-worker restart fixture passed with this wrapper: 13 genes,
three families, 42 directed ortholog pairs and one expected duplication.
Every known family produced exactly one start and one return trace. The
missing-ID fixture still exits zero without complete outputs and is detected
by the independent fixture auditor. Fixture receipts are archived as
`metadata/orthology_native_traced_restart_fixture_receipt.json` and
`metadata/orthology_native_traced_restart_fixture_readback.json`.

Recovery uses fresh physical input copies under
`results/orthology/expanded-reconciliation-execution-v2/`, with both complete
guide partitions preserved. The plan is
`metadata/expanded_reconciliation_execution_recovery_plan.json` and service
`fungal-expanded-reconciliation-recovery-20260922.service`.
The existing eight-CPU quota, four reconciliation workers, 128 GiB memory
limit, no swap, file-descriptor allowance, pre-run memory/disk gates and
1 TiB emergency disk guard remain in force. Existing storage/runtime
estimates are retained; the longer stall allowance is not a completion ETA.
No GPU prediction or paid resource is used. Full output validation and the
second guide execution remain outstanding.
