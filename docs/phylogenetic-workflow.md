# Sequence-based phylogenetic framework

The first backbone uses the 125 eukaryota_odb12.2 BUSCO markers, allowing the fungal ingroup and all 25 outgroups to enter the same marker screen. This is a starting marker panel, not evidence that the eventual tree is robust. Broad single-copy status does not exclude ancient paralogy, compositional bias, contamination, or model misspecification.

After successful BUSCO QC for the complete analysis manifest:

```bash
python scripts/prepare_phylogenetic_markers.py --output results/phylogeny/markers-full-v1
python scripts/align_phylogenetic_markers.py --markers results/phylogeny/markers-full-v1 --output results/phylogeny/alignments-full-v1
```

Extraction requires successful BUSCO receipts tied to unchanged input FASTAs. It retains only Complete single-copy calls, validates marker membership and exact protein sequence identity, and records protein IDs, checksums and occupancy. Duplicated and fragmented hits are not resolved by arbitrarily taking one sequence. Each taxon is represented once per retained marker. Missingness is measured against all planned taxa, with ingroup and outgroup counts reported separately.

The optional `--allow-incomplete` extraction switch supports inspection of completed jobs during the full-dataset run. Its receipt identifies the snapshot as incomplete, lists pending taxa, and prevents downstream alignment. It is not a pilot dataset or a smaller substitute for the full analysis. Extraction output directories are immutable; use a new name when adding taxa.

Alignment uses MAFFT 7.525 `--auto --inputorder`, four concurrent alignments with two threads each by default. This is eight CPU threads on the existing host, no GPU or paid resources. The initial 125-marker panel should require modest storage relative to the proteomes; runtime will be measured during the full run because divergent, long markers may be much slower. Record MAFFT's selected algorithm in each log. Output validation requires unchanged taxon identities and ungapped residues. Resume verifies input/output checksums, command and version. Very low occupancy markers with fewer than four sequences are explicitly skipped. All other occupancy filtering remains an analytical choice to be documented before tree inference.

## Required steps after alignment

1. Quantify occupancy by lineage, sequence length outliers, alignment uncertainty, composition and unusually long branches. Evaluate annotation or contamination explanations for outliers before removal, preserving exclusions and reasons.
2. Define a primary filtering policy and retain unfiltered and alternative masks for sensitivity analyses. Avoid filtering solely to improve support for a preferred topology.
3. Infer individual gene trees and a partitioned concatenated backbone with support estimates. Compare a gene-tree species estimate, account for gene-tree estimation error, and summarize concordance/conflict. Deep splits require models and sensitivity checks appropriate to compositional and site heterogeneity; a single homogeneous model is insufficient evidence.
4. Add a richer fungal marker panel for the ingroup. Retain shared markers and suitable close outgroups to connect this analysis to the broad backbone. Compare marker panels, reduced distant-outgroup sets, long-branch exclusions and the narrower fungal circumscription.
5. Establish the species tree's uncertainty before gene-tree reconciliation and branch-specific sequence/structure models. Do not use structural clusters as orthogroups or treat structural-alphabet distances as additive physical displacements.

These downstream inference steps have not yet been executed. Time calibration, orthogroup gene trees and reconciliation remain separate deliverables.
