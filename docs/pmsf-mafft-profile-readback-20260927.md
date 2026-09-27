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
