# Phylogeny interruption recovery, September 28

The 12:28 EDT interruption also affected the fourth crossed PMSF run and
the MAFFT marker-tree workflow. The PMSF controller journal records that its
runner died with SIGTERM. The MAFFT producer and audit PIDs were absent;
the downstream topology controller failed because its required full-audit
receipt did not exist. The source of the termination is unresolved.

## Marker trees resumed

Before restarting, the original plan's input/software hashes were checked,
51 completed marker receipts and tree hashes were verified, all available
gzip checkpoints were read, and 702 files were copied into
`results/recovery-20260928/mafft-markers-before-resume` with a checksum receipt.
No replacement marker or PMSF native processes were found during inspection.
The original marker controller was restarted with the unchanged plan,
four two-thread workers, 24 GiB memory limit and no swap. IQ-TREE's existing
checkpoints are retained; completed marker results are reused by the runner.
These are execution-level checks, not the complete input/support audit.

New launch identities and plans reconnect the full marker audit and subsequent
alignment-topology comparison. Previous failed controller directories are
preserved; new controller outputs use September 28 paths. All three resumed
controller identities were verified live after launch.

## Fourth PMSF run prepared, not restarted

`metadata/pmsf_checkpoint_recovery_plan_20260928.json` binds the unchanged
input/software configuration and readable native checkpoint. The checkpoint,
configuration and logs are separately preserved under
`results/recovery-20260928/pmsf-before-resume`.
The original controller/native PIDs were absent and no matching PMSF process
was active. The log ended during normal NNI optimization, without a completion
receipt. The recovery is for a process interruption on the same boot.

At preparation about 690 GiB memory was available, below the existing 750 GiB
launch prerequisite, so the recovery was not launched. Recheck live process
identities, checkpoint/input hashes and memory before executing the prepared
plan with `scripts/recover_species_pmsf_checkpoint.py`. Preserve the 16-thread,
650 GiB service-memory, no-swap limits. After launch, reconnect the full PMSF
readback to the new controller identity using a new plan/output. Do not restart
the old crossed controller directly: it rejects the unfinished output directory.

The project remains incomplete. Neither missing receipts nor inactive old
processes imply successful scientific results.

## Memory-gated recovery queued

The v2 recovery plan and `scripts/wait_for_pmsf_recovery_memory.py` now queue
the restart automatically. The wrapper checks available memory every 30 seconds,
verifies the pinned recovery dependencies, and refuses a duplicate native run
or an already produced receipt. The existing recovery runner then repeats its
checkpoint/input, memory, disk and global-lock checks immediately before launch.
The 750 GiB prerequisite is unchanged. The wrapper was verified live and waiting
at 687.20 GiB available; it has not started native inference at that observation.

`fungal-pmsf-memory-recovery-20260928.service` retains the 16 CPU, 650 GiB memory
and no-swap limits. Its new process identity is pinned in
`metadata/pmsf_memory_recovery_launch_20260928.json`. The full profile/tree/1,000
bootstrap readback is queued against that identity with the new
`metadata/pmsf_fourth_readback_plan_20260928_v2.json`. The old recovery and
failed readback records remain intact. No GPU jobs are launched by this queue.
