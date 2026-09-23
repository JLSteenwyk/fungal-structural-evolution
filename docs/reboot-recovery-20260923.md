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
