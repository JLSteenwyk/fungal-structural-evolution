# September 23 checkpoint recovery

The host booted at 09:20 Eastern on September 23. The previous journal ends
at September 22 15:58:28; that is the last recorded message, not proof of a
specific shutdown cause. Previous transient project services are absent after
reboot and the original job processes are gone. Saved controller states that
say `running` are stale and are not evidence of active computation.

The previous turn was interrupted while examining root HOG coverage. No root
completeness result was produced. In particular, the source contains 500,113
singleton genes, 200,703 genes in two-/three-gene families and 5,115,031 genes
in larger families. The 5,274,605 observed N0 assignments do not cover all
5,315,734 nonsingleton genes; the 41,129-gene difference needs examination of
native HOG semantics and gene disposition before any completeness claim.

## Restored and verified live

- The third crossed species-tree run (`pmsf-mafft-profile-v1`) resumed with
  the original IQ-TREE executable, inputs, settings and output prefix.
  Its compressed checkpoint passed integrity reading and was copied to
  `results/recovery-20260923/pmsf-mafft-profile/original_checkpoint.ckp.gz`
  before native execution. IQ-TREE explicitly restored the site-frequency
  model, initial tree and 1,000 bootstrap trees. The 16-thread/600G native
  request remains within the original resource plan; the service has a
  650 GiB memory limit and no swap. Completion will still require full
  output/profile/support verification.
- The 122-marker paired fitting run resumed unchanged. All four interrupted
  AA fits emitted checkpoint-restoration messages. Four single-thread workers
  retain the original 2 GiB per-fit request; the service limit is 16 GiB.
  Source/model/configuration bindings are checked, and a full fit audit runs
  after inference. Recovery preserves the original plan and adds explicit
  recovery provenance to its eventual completion receipt.
- The full profile-guide HOG identity audit restarted read-only into fresh
  `results/orthology/profile-hog-identities-20260923-v3`. The interrupted v2
  stopped after logging N108 and has no final receipt. It is preserved.
  The new one-CPU audit checks all 525 tables from the beginning.
- The fourth crossed species-tree run is queued behind the new exact species
  recovery process identity. The 73,200-draw paired-resampling stage is queued
  behind the new exact fit-recovery identity. Both still require successful
  bound completion receipts before starting native work.

The new recovery scripts are `scripts/recover_species_pmsf_checkpoint.py` and
`scripts/recover_paired_fits_checkpoint.py`. Plans and observed live process
identities, child commands, resource limits and boot identity are recorded in
`metadata/recovery_20260923_*`. The PMSF runner retains the existing global
species-tree lock; the paired fitting runner retains its existing output lock.
No force-restart flag, changed scientific model or GPU inference was used.
Recovery runtime is not a new measured completion ETA. Original plans remain
in force; the work uses existing authorized hardware without new charges.

## Recovery still required

The MAFFT-guide reconciliation and its both-guide tree checks/small-family
supplement, whole-proteome Foldseek database construction, refreshed AlphaFold
residue/features/confidence/paired-fit chain, PAE prefetch and background
retrieval remain stopped pending checkpoint/output review and new process
bindings. Their old PIDs must not be used as live dependencies. The completed
catalogs, verified ESMFold inputs, first-guide native output and archived
receipts are retained. The unrelated IQ-TREE process already running after
boot was left untouched. GPU structure prediction remains paused.

## Structural atlas and AlphaFold chain restored

A second recovery pass has verified seven live services. The original completed
catalog and mapping are reused with their existing checksums. Every interrupted
output directory is preserved; replacement runs use September 23 output paths.
No checkpoint-resume capability is assumed for incomplete Foldseek database
construction or residue-level readback: those stages run again from their
verified inputs.

| Stage | CPU quota | Memory limit | Recovery status |
|---|---:|---:|---|
| Whole-proteome Foldseek database | 4 | 64 GiB | Running, fresh database output |
| Refreshed AFDB residue readback | 1 | 16 GiB | Running, complete mapping reused |
| AFDB PAE prefetch | 2 | 8 GiB | Running, validated global cache reused |
| Native structural features | 4 | 16 GiB | Waiting for new residue auditor |
| Mapping-bound confidence and full-context readback | 2 | 16 GiB | Waiting for new feature/prefetch processes |
| Paired AA/3Di input preparation/readback | 1 | 32 GiB | Waiting for new confidence process |
| Supported paired point fits and audit | 4 | 16 GiB | Waiting for new paired-input process |

All services have no swap and reduced CPU/IO priority. Existing per-stage disk,
RAM and runtime planning limits are retained in their recovery plans; quotas are
limits rather than reserved CPU usage while waiting. The new process identities,
creation times/start ticks, command lines and second live-status observations
are recorded in `metadata/recovery_20260923_afdb_launch_inventory.json`.
Dependency plans bind the new producer receipts and plan checksums. The completed
mapping's original preparation plan remains unchanged so its receipt retains
its exact original provenance.

PAE prefetch uses the existing `data/structures/pae` cache. Cached receipts and
matrix bytes must pass the existing checks before reuse. A replacement manifest
will be produced rather than accepting the interrupted manifest. Foldseek output
is now `results/structural_clusters/whole-proteome-afdb-database-20260923-v1`;
the fixed source catalog still contains 1,290,278 models. No structural search,
clustering or evolutionary result is claimed by merely restoring these jobs.

This supersedes the stopped status for these stages in the initial recovery
inventory above. MAFFT-guide reconciliation and its dependent checks/small-family
supplement, plus the background full-catalog retrieval job, still need recovery.
The species-tree, marker fits, HOG verification and their reconnected downstream
jobs remain separate live recovery branches. GPU prediction remains paused.

## Reconciliation and background retrieval restored

The remaining reconciliation and retrieval recovery services are now live.
`scripts/recover_expanded_reconciliation.py` preserves the original controller's
native execution and verification logic while permitting explicit reuse of a
completed guide. The profile stage's archived completion status, original plan
binding, mandatory output hashes and every copied source-file hash are checked
before reuse. A new directory link exposes that preserved profile result to
read-only downstream checks; no profile native inference is rerun. The MAFFT
stage receives fresh physical input copies and fresh native output under
`results/orthology/expanded-reconciliation-execution-20260923-v3`.

The incomplete previous MAFFT output remains intact. No resumable native
reconciliation checkpoint was established, so the replacement reruns the stage.
The same eight CPU equivalents, four native analysis workers, 128 GiB memory,
no swap, six-hour task-completion allowance and disk/memory gates are retained.
Five completion-gate checks reject incorrect guide, incomplete status, changed
source flag, missing mandatory artifacts and incorrect original-plan binding.
Both-guide resolved-tree validation and the small-family orthologue supplement
are queued against the replacement controller's exact PID/start time and plan.
All mandatory outputs are rehashed before the eventual combined completion
handoff. Full semantic reconciliation validation remains outstanding.

The unchanged background `retrieve_matched_models.py` also resumed through
`scripts/recover_structure_retrieval.py`. Before launch, every record of the
1,557,799,515-byte JSONL log was parsed: 1,329,102 records, 1,326,490 distinct
latest accessions, including 1,322,928 latest verified statuses. These are log
statuses, not newly audited unique models. The original producer rehashes
cached coordinate bytes before constructing its remaining queue, then retains
the existing two HTTP-worker limit and exact full-sequence checks. It keeps
its exclusive cache lock, appends new records and updates the unversioned live
queue snapshot. Frozen catalogs used by running analyses remain unchanged.

Retrieval recovery is capped at two CPU equivalents and 64 GiB RAM with no
swap and reduced scheduling priority. It starts only with at least 3 TiB free;
a wrapper requests a clean stop below 1 TiB, allowing outstanding requests to
drain before escalation if needed. The 2 TiB output and 24–1,008-hour runtime
allowances are conservative, uncalibrated planning bounds for the remaining
full catalog, not a new completion ETA or a commitment to fill every missing
protein. No new paid infrastructure or GPU prediction is involved.

The replacement plans, initial launch records and a second live process check
are in `metadata/recovery_20260923_reconciliation_*`,
`metadata/recovery_20260923_{small_family_supplement,resolved_tree_readback}_plan.json`,
`metadata/recovery_20260923_structure_retrieval_*` and
`metadata/recovery_20260923_remaining_jobs_verification.json`.
This supersedes the pending recovery status for these jobs above. Restoration
of execution is not completion of the scientific analyses.


## PAE cache repair and dependent-controller replacement

The refreshed PAE prefetch terminated with 30,587 verified matrices and one
failure. Both cached files for AF-A0A0C3AZU7-F1 version 6 (compressed matrix
and receipt) were zero bytes. They were preserved under
`results/recovery-20260923/pae-empty-cache-quarantine/` while holding the cache
lock. A fresh version-matched retrieval passed matrix validation (340 residues).
The cause of the empty files has not been established.

`scripts/recover_refreshed_pae_manifest.py`, with
`metadata/refreshed_afdb_pae_repair_plan.json`, rebuilt a fresh complete manifest
at `results/structural_pae/gdm-current-prefetch-repaired-20260923-v1`.
Every compressed matrix was rehashed; all 30,587 previously successful records
were unchanged, and model/version/sequence/length/URL bindings were checked
against the full catalog. All 30,588 matrices are now verified, with zero
failures. Numerical validation was reused for unchanged matrices; complete
mapping-bound validation and directional context readback remain downstream.
The completed result and controller receipts are archived as
`metadata/refreshed_afdb_pae_repair_{completed,controller}_receipt.json`.

Three waiting controllers were confirmed live with no children, stopped, and
replaced: confidence qualification, paired input preparation, and paired fits.
The new plans are `metadata/pae_repair_20260923_afdb_*_plan.json`; process
identities and exact launch commands are in
`metadata/pae_repair_20260923_afdb_launch.json`. All new dependent outputs use
fresh v2 paths. The healthy native coordinate audit and its v1 output are
preserved. The completed PAE repair is bound by plan and result checksums;
PID zero explicitly represents an already completed predecessor.

Live replacement process identities and all cross-stage paths were checked
after launch. The confidence controller waits for the ongoing coordinate
audit; the other two controllers wait on their exact new predecessors.
Resource limits remain 2 CPU/16 GiB, 1 CPU/32 GiB, and 4 CPU/16 GiB respectively,
with no swap or GPU use. This repairs the pipeline dependency; it does not
establish completion of confidence qualification or evolutionary fitting.
