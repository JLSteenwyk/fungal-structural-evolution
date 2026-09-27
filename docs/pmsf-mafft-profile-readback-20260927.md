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
