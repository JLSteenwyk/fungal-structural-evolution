# Complete marker-tree diagnostics, 22 September 2026

All 125 planned marker gene trees are inferred and audited. This supersedes
the earlier 70-marker snapshot for full-batch diagnostics; the earlier outputs
remain preserved. Independent reconstruction checked 23,441,199 retained
alignment characters and 59,315 internal support splits. Support labels are
SH-aLRT from 1,000 replicates, not bootstrap percentages.

The complete guide comparison covers 261,500 marker/guide/edge/cutoff rows:
125 markers, two homogeneous guides, 523 internal guide edges and two descriptive
SH-aLRT cutoffs (80 and 95). Independent readback verified the complete grid,
restricted splits, support labels, classifications, conflict witnesses and
2,092 edge summaries. Exhaustive independent searches for maximum conflicting
support covered five deterministic guide edges per marker/guide (1,250 cases);
other maxima were not independently exhaustively searched.

Both homogeneous guides contain the exact split separating the 501 fungal
entries from the 25 outgroups. Among marker trees, 21 contain the role-separating
split and 104 do not. Of those 21, 17 have SH-aLRT support at least 80 and 11 at
least 95. Each marker is assessed on its own retained taxon universe, requiring
at least two entries from each role; all 125 are informative by this criterion.
Independent graph traversal checked both guides and all marker trees.

Absence of the exact split does not establish a biological cause or prove
confident non-monophyly. Gene-tree error, missing taxa, model inadequacy,
annotation/copy issues and biological discordance require further evaluation.
An unrooted separating edge does not identify a root coordinate or establish
event direction. Entry-level taxonomy uncertainty remains. The guide trees
are provisional homogeneous alternatives, not the final supported species tree.
Restricted guide splits may represent collapsed paths; rows are not independent
branch replicates, gene concordance factors or tests of structural acceleration.

These diagnostics provide the complete marker set for the next phylogenetic
sensitivity and branch-eligibility analyses. Supported mixture-tree comparisons,
rooting sensitivity, reconciliation validation and uncertainty-aware structural
branch tests remain unfinished.

## Reproducibility

The five-stage controller is `scripts/run_complete_marker_diagnostics.py`.
`metadata/complete_marker_diagnostics_plan.json` pins its input files, scripts,
commands, expected completion states and resources (one CPU, 16 GiB memory,
no swap, 10 GiB output allowance, 50 GiB free-disk gate; 0.5–8 hours was an
uncalibrated planning range). Execution completed without a GPU.

The controller receipt is `metadata/complete_marker_diagnostics_receipt.json`;
all five stage receipts are archived as `metadata/complete_marker_*_receipt.json`.
The small per-tree role-separation table is
`metadata/complete_marker_root_split_diagnostics.tsv`. Large support and conflict
tables remain at the result paths recorded in the plan and receipts. All five
stage receipt hashes and their declared artifact hashes were checked before
archiving. Reproduction requires fresh output paths and a newly pinned plan;
the controller refuses to overwrite existing results.
