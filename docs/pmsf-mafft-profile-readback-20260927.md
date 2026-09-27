# MAFFT-alignment PMSF run with profile-derived guide

The third crossed species-tree run completed and passed the full saved-output
readback on September 27. It uses the 63,750-site MAFFT alignment, all 526 taxa,
and the profile-alignment guide. The fourth combination (MAFFT alignment and
MAFFT-derived guide) has started through the existing serial controller.

All 63,750 ordered amino-acid frequency vectors were checked for dimensions,
finite bounded values and sums at output precision. All 1,000 bootstrap trees
were checked for exact tips and finite nonnegative branches. Their empirical
split frequencies agree with the NEXUS weights and printed tree support at
reported precision. The bootstrap ensemble contains 1,230 distinct splits,
including terminal splits. Among 523 internal maximum-likelihood branches,
464 meet the previously used descriptive label of SH-aLRT >=80 and empirical
UFBoot >=95.

The maximum-likelihood and reoptimized consensus trees have RF distance four:
two internal splits in each tree differ. Both trees and explicit split support
remain available; neither is selected here as the final framework. The reported
consensus likelihood improvement is 0.267054. Likelihoods were read back, not
independently recomputed, and SH-aLRT was not rerun. The report repeats a warning
that 31 sequences contain more than 50% gaps or ambiguity; three warning lines
are not three independent groups of 31 taxa.

Reproduce the audit with a new output directory:

```bash
python scripts/audit_species_pmsf.py --run results/phylogeny/pmsf-mafft-profile-v1 --output NEW_AUDIT_DIRECTORY
```

Completed audit: `results/phylogeny/pmsf-mafft-profile-readback-20260927-v1`.
Completion proof: `metadata/pmsf_mafft_profile_readback_completed_20260927.json`.
The complete inference and audit artifact hashes were rechecked after successful
audit termination. The audit consumed about 16.6 CPU seconds and peaked at
98.3 MB within its one-CPU/8-GiB allowance. Inference was not rerun.

The fourth crossed run, cross-alignment comparisons, taxon/marker sensitivity,
rooting, model adequacy and discordance remain necessary. Neither a support
threshold nor a successful output audit establishes biological correctness.

## Three-run topology comparison

All three validated ML trees were compared by canonical unrooted splits. An
independent DendroPy parser and bipartition-compatibility check reconstructed
all 1,662 split-presence rows and 132 incompatible split-pair records. The three
trees share 493 of their 523 internal splits.

| Runs (alignment / guide) | Shared internal splits | RF distance | Incompatible pairs meeting both support criteria |
|---|---:|---:|---:|
| Profile/profile vs profile/MAFFT | 516 | 14 | 0 |
| Profile/profile vs MAFFT/profile | 493 | 60 | 1 |
| Profile/MAFFT vs MAFFT/profile | 499 | 48 | 1 |

The two supported-conflict records represent the same distinct split conflict,
not two independent biological events. Both profile-alignment trees group
Saccharomyces arboricola (F706196) with the six shared entries F114524, F114525,
F27291, F27292, F332112 and F4932. The MAFFT-alignment tree instead places
S. eubayanus (F1080349) and S. uvarum (F230603) with those six in the conflicting
split. The shared entries include S. pastorianus and the named S. cerevisiae ×
S. kudriavzevii hybrid. This overlap with hybrid entries motivates the already
required hybrid-exclusion sensitivity; it does not establish hybridity as the
cause of the alignment-dependent conflict.

For the conflicting profile split, SH-aLRT/empirical UFB support is 97.6/99.8
with the profile guide and 89.7/99.5 with the MAFFT guide. The conflicting
MAFFT/profile split has support 99.3/100.0. Every conflict, including those below
the existing SH-aLRT80/UFB95 labels, remains in the complete results. Raw
likelihoods and site-frequency vectors are not compared across alignments
because the site universes differ. These are unrooted split comparisons, not
claims about rooted clades or evolutionary direction.

Reproduction scripts: `compare_completed_species_pmsf.py` and
`readback_completed_species_pmsf_comparison.py`. Full tables are under
`results/phylogeny/pmsf-three-run-topology-sensitivity-20260927-v1`; receipts,
independent proof and compact comparisons are archived in
`metadata/pmsf_three_run_topology_*_20260927.*`. The fourth run and consensus,
rooting, taxon/marker and model-adequacy sensitivities remain unfinished.

## Marker-level assessment of the Saccharomyces conflict

All 125 audited profile-alignment marker trees were compared with both
conflicting full species-tree splits, restricted to each marker's observed
taxa. The saved trees were also projected after removing hybrid tips F27292
and F332112. This is not a fresh hybrid-excluded alignment, inference or support
analysis. Where tip removal collapses several edges to one split, the diagnostic
uses the minimum source SH-aLRT only when every source support is available.
Missing support remains unresolved.

| Saved-tree condition | Alternative split | SH-aLRT cutoff | Concordant | Conflicting | Unresolved | Uninformative coverage |
|---|---|---:|---:|---:|---:|---:|
| All available taxa | Profile | 80 | 18 | 58 | 46 | 3 |
| All available taxa | MAFFT | 80 | 24 | 59 | 39 | 3 |
| All available taxa | Profile | 95 | 5 | 14 | 103 | 3 |
| All available taxa | MAFFT | 95 | 19 | 15 | 88 | 3 |
| Hybrid tips removed | Profile | 80 | 17 | 53 | 52 | 3 |
| Hybrid tips removed | MAFFT | 80 | 23 | 59 | 40 | 3 |
| Hybrid tips removed | Profile | 95 | 4 | 10 | 108 | 3 |
| Hybrid tips removed | MAFFT | 95 | 18 | 13 | 91 | 3 |

These counts demonstrate disagreement and limited resolution among the saved
marker trees. Conflict with one candidate need not imply support for the other;
other arrangements are possible. Counts are not independent votes, gene
concordance factors, probabilities or a causal explanation. The same markers,
including their missing-data and copy caveats, recur across conditions. All
marker trees here use the profile alignment, so this does not replace a matched
MAFFT-marker comparison or fresh hybrid-excluded inference.

The complete 1,000 assessment rows and eight summaries were independently
reconstructed directly from gene-tree files using DendroPy and bitmask split
compatibility. Scripts are `assess_saccharomyces_marker_conflict.py` and
`readback_saccharomyces_marker_conflict.py`. Full results are in
`results/phylogeny/saccharomyces-marker-conflict-20260927-v1`; metadata archives
include receipt, independent proof and summary under
`saccharomyces_marker_conflict_*_20260927.*`.

## Fresh full-matrix hybrid-excluded inference started

Prepared both complete species matrices with only F27292 and F332112 removed:
524 taxa remain, including all 25 outgroups. The profile and MAFFT matrices
retain all 49,027 and 63,750 original columns, respectively. An independent
FASTA parser checked every one of 59,095,148 retained characters and the exact
taxon difference; partition and site-mapping files are byte-identical to their
sources. Input preparation is reproducible with
`scripts/prepare_hybrid_excluded_species_matrices.py`; receipt and readback are
archived under `metadata/hybrid_excluded_matrix_*_20260927.json`.

Fresh full-matrix guide inference is now running under LG+F+G4, the same
homogeneous model used for the original guides. Both matrices run serially with
four threads, IQ-TREE memory 32 GiB, service memory 40 GiB and no swap. Each has
1,000 SH-aLRT and 1,000 UFB replicates with bootstrap NNI and saved bootstrap
trees. Identical sequences are retained. The seed is 20260927. No original
matrix or inference is overwritten. The controller can continue an identical
interrupted stage using IQ-TREE's checkpoint; it checks all pinned inputs before
and after execution.

Resource planning reserves 10 GiB output and requires 100 GiB free disk and
64 GiB available RAM before each stage. The broad allowance is 4–168 hours per
matrix, not a measured ETA. For context, original 16-thread guides without
support took approximately 15.6–15.8 wall hours and 226–230 CPU hours; these
four-thread supported runs can take substantially longer. They use existing
CPU resources and incur no new charges. The existing fourth full-taxon PMSF
run continues separately.

Plan and live identity: `metadata/hybrid_excluded_species_guide_plan_20260927.json`
and `metadata/hybrid_excluded_species_guide_launch_20260927.json`.
Runner: `scripts/run_hybrid_excluded_species_guides.py`.
Output: `results/phylogeny/hybrid-excluded-guides-20260927-v1`.
These are full-sampling homogeneous-model sensitivities and fresh guide inputs,
not completed results or substitutes for subsequent mixture-model analyses.
Full result/support audits, mixture-model guide sensitivity, and the other
identity, contamination and marker-sampling checks remain required.

## Fourth-run full readback queued

The fourth (MAFFT/MAFFT) run now has an automatic full saved-output audit
queued behind its exact live controller identity. The new controller
`scripts/advance_fourth_pmsf_readback.py` checks the producer PID, creation time
and command, waits for termination, requires a successful completed inference
receipt, then runs the unchanged `audit_species_pmsf.py`. Inputs, execution
configuration and audit scripts are pinned. The audit requires all 526 taxa,
63,750 profiles and 1,000 bootstrap trees, and rechecks output hashes.

Resource allowance is one CPU, 8 GiB RAM, no swap and 100 MB output; the
0.01–2-hour post-inference planning range is conservative relative to the third
run's 16.6-CPU-second readback. This starts no additional inference. Plan and
verified live identity are archived as
`metadata/pmsf_fourth_readback_{plan,launch}_20260927.json`. Actual results are
pending under `results/phylogeny/pmsf-mafft-mafft-readback-20260927-v1`, with a
separate handoff receipt directory. All-four comparison and remaining framework
sensitivities still require completion after the audit.

## Full MAFFT marker-tree counterpart started

All 125 marker trees are now being inferred from the complete MAFFT species
matrix. The existing profile-marker coverage rule is retained: at least 50
unambiguous amino acids and at least 30% of a marker's columns per taxon.
The full planning grid independently reconstructed from matrix characters has
418–504 retained taxa per marker and covers all 63,750 columns. This is the
complete marker batch, not a selected Saccharomyces subset. Both hybrid entries
remain available here; this is alignment-method sensitivity, separate from the
hybrid-excluded species runs.

A separate runner preserves the original profile jobs and uses the same
marker-hash seeds, model search (LG/WAG/JTT with +F and gamma), and 1,000 SH-aLRT
replicates. Four concurrent jobs each use two CPU threads and at most 4 GiB
IQ-TREE memory. The service has 24 GiB RAM, no swap, and checks 64 GiB available
RAM and 100 GiB free disk before beginning. Output allowance is 10 GiB. The
12–168-hour batch planning envelope is not a guaranteed ETA; original profile
marker fits summed to about 250 case-hours with a 76-minute median, and MAFFT
markers have more columns. Existing CPU resources incur no new charges.

Runner: `scripts/run_mafft_marker_gene_trees.py`.
Output: `results/phylogeny/mafft-marker-gene-trees-20260927-v1`.
Plan, independent planning-grid readback and verified controller/four-worker
identities are archived under `metadata/mafft_marker_tree_*_20260927.*`.
Actual emitted input reconstruction, full support readback, copy/coverage
review and common-taxon comparisons remain required. SH-aLRT values are not
bootstrap proportions or posterior probabilities, and no gene-tree result is
claimed complete at launch.

The full MAFFT marker batch now has an automatic two-stage audit queued behind
its identified producer. First, a separate matrix-parameterized copy of the
existing support auditor records every internal split and its reported support.
Second, the unchanged independent input/graph checker reconstructs every marker
alignment from the full source matrix and verifies the split identities by
undirected graph-edge removal. The original profile-audit script remains intact.
The controller requires all 125 successful fits and rejects incomplete batches.

Scripts: `audit_mafft_marker_tree_support.py` and
`advance_mafft_marker_tree_audit.py`; independent checker:
`readback_marker_tree_snapshot.py`. The pinned plan and verified waiting-process
identity are archived in `metadata/mafft_marker_tree_audit_{plan,launch}_20260927.json`.
The audit uses one CPU, 8 GiB RAM, no swap and a 1-GiB output allowance, with a
0.1–4-hour planning envelope after inference. Outputs will be
`results/phylogeny/mafft-marker-support-complete-20260927-v1` and
`results/phylogeny/mafft-marker-support-readback-20260927-v1`. No pending output
is claimed verified; biological adequacy and common-taxon comparisons remain
separate downstream requirements.

## Exact marker coverage comparison across alignments

All 65,750 marker–taxon cells (125 markers × 526 taxa) were reconstructed from
both full matrices under the same coverage rule. Both alignments retain 59,690
cells; MAFFT alone retains 149; profile alone retains zero; neither retains
5,911. Taxon membership differs for 30 markers. Thus all profile-marker input
taxa are contained in the planned MAFFT input set for the corresponding marker.
Common sets contain 418–504 taxa, and exact common/additional lists are retained
for every marker. The 149 additions are marker–taxon combinations, not 149
newly sampled species.

`scripts/compare_marker_alignment_coverage.py` constructs the complete ledger
and checks profile membership against finished inference inputs and MAFFT
membership against its complete planning grid. A separate FASTA parser and
NumPy residue counter in `scripts/readback_marker_alignment_coverage.py`
reconstruct every coverage cell and exact membership list. Full results are
in `results/phylogeny/marker-alignment-coverage-20260927-v1`; summary, receipt
and independent proof are archived under
`metadata/marker_alignment_coverage_*_20260927.*`.

These sets prepare common-taxon topology comparisons after the MAFFT batch
passes its output audit. Actual MAFFT emitted inputs and tree tips still must
match the planned sets. Pruning outputs to shared taxa will control the
comparison's taxon universe, but cannot undo how additional taxa influenced
inference. A matched-input refit would address that separate sensitivity.
