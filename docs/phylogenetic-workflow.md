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

## Alignment audit and profile comparison

Run `python scripts/assess_marker_alignments.py --alignments results/phylogeny/alignments-full-v1 --output results/phylogeny/alignment-audit-full-v1` after the MAFFT batch finishes. Inspection of a live batch is supported with `--allow-incomplete`, which labels the output staging. Audits report unambiguous-residue occupancy, constant and parsimony-informative columns, per-taxon retained residues, and explicit 1-based column masks at 25%, 50% and 75% occupancy among taxa present for each marker. Ambiguous residues do not count as observed amino-acid states. These masks quantify missingness; they are not probabilistic alignment-confidence estimates or a final filtering decision.

An executed audit of the first 11 completed MAFFT alignments found median retained fractions of 8.6%, 7.4% and 6.8% under these masks. Full protein lengths can include lineage-specific extensions: for marker 115730at2759 the median protein length is 358 aa, maximum 2,058 aa, and MAFFT alignment length 4,061 columns. Large insertion-rich regions motivate a second full-panel alignment approach; they are not by themselves grounds for discarding taxa.

```bash
python scripts/align_markers_to_profiles.py --markers results/phylogeny/markers-full-v1 --output results/phylogeny/profile-alignments-full-v1
python scripts/assess_marker_alignments.py --alignments results/phylogeny/profile-alignments-full-v1 --output results/phylogeny/profile-alignment-audit-full-v1
```

The profile comparison uses HMMER 3.4 `hmmalign` against the pinned BUSCO HMMs with four concurrent processes on the existing host. It preserves full Stockholm alignments, verifies ungapped sequence identity, and extracts profile match-state columns according to the reference annotation, requiring their count to equal the HMM length. Profile, input, full-alignment and extracted-alignment checksums are retained, together with the original alignment column indices. Stockholm posterior-probability annotations are preserved for later assessment. This profile-based approach also has model assumptions and does not establish a uniquely correct alignment; compare phylogenetic sensitivity against MAFFT and inspect leading conflicts.

## Initial concatenated guide

```bash
python scripts/build_species_matrix.py --alignments results/phylogeny/profile-alignments-full-v1 --audit results/phylogeny/profile-alignment-audit-full-v1 --output results/phylogeny/profile-matrix-50-v1
python scripts/run_initial_species_tree.py --matrix results/phylogeny/profile-matrix-50-v1 --output results/phylogeny/initial-guide-v1
```

The initial matrix uses profile columns with at least 50% unambiguous occupancy among taxa present per marker. This is a topology-blind baseline; alternative occupancy masks and MAFFT-based alignments remain required. Concatenation recomputes masks, inserts gaps for absent taxa, normalizes nonstandard/ambiguous residues to X, and preserves a site-to-marker coordinate table and a NEXUS partition definition. It rejects taxa with no unambiguous residues. The executed matrix contains 526 taxa and 49,027 columns.

The launched IQ-TREE 3.0.1 guide uses unpartitioned LG+F+G4, seed 20260913, 16 threads and a 32 GB limit. It omits bootstrap support and does not establish the final species tree. IQ-TREE estimated approximately 15.9 GB memory and reported nominal composition-test failures for 520/526 sequences; this reinforces the need for richer models and compositional sensitivity rather than interpreting the homogeneous-model guide as robust. No taxa were removed based on these tests.

Coverage needs explicit sensitivity analysis: Amoeboaphelidium protococcarum contributes only 164 unambiguous positions because its raw-proteome BUSCO results were 95.2% duplicated and 0.8% single-copy. Its annotation maps every protein to a distinct gene feature, so ordinary isoform collapsing does not resolve this issue. Pirum and Abeoforma contribute 1,224 and 1,804 positions, respectively. Their placements cannot be treated as well-established solely because a tree contains those tips. Gene-tree/paralogy and assembly-quality investigation remain necessary.

## Individual marker trees

Run `python scripts/run_marker_gene_trees.py` on the frozen profile-matrix-50-v1 inputs. It processes all 125 markers, retaining a taxon within a marker only when it has at least 50 unambiguous amino acids and at least 30% of that marker's retained columns. Exclusions and observed residue counts are recorded for every marker; this does not remove taxa from the project or species matrix. The rule is a coverage baseline requiring sensitivity assessment, not evidence that shorter sequences lack phylogenetic information.

Four concurrent IQ-TREE jobs, two threads each, compare LG/WAG/JTT with empirical frequencies and gamma-distributed rates. ModelFinder's selected model and criteria are retained in each IQ-TREE report. Trees receive 1,000 SH-aLRT replicates; these are approximate likelihood-ratio supports, not bootstrap supports. Fixed marker-derived seeds and input/configuration/tree checksums support reproducibility. IQ-TREE checkpoint resume is allowed only under unchanged configuration. Pending jobs are cancelled if a job reports failure.

The live run is `results/phylogeny/marker-gene-trees-v2`. Version v1 is a terminal failed launch: IQ-TREE rejected `--mfreq`; v2 uses the documented single-dash `-mfreq` and `-mrate G` spellings and has entered likelihood optimization. Preserve v1 logs as failure evidence. Completed trees must still be checked for gene-tree estimation error, support-sensitive discordance and compositional/model sensitivity before reconciliation or species-tree conclusions.

## Full MAFFT and profile correspondence sensitivity

All 125 full-protein MAFFT alignments completed and passed residue/identity preservation checks. The complete audit contains 697,866 raw columns, of which 63,750 pass 50% occupancy among taxa present in each marker. The profile alignments contain 50,941 HMM match-state columns, of which 49,027 pass the same occupancy rule. These different alignment scopes preclude interpreting raw retention fractions as relative accuracy.

Reproduce the completed audit and alternative matrix using new output directories:

```bash
python scripts/assess_marker_alignments.py --alignments results/phylogeny/alignments-full-v1 --output results/phylogeny/mafft-audit-full-v1
python scripts/build_species_matrix.py --alignments results/phylogeny/alignments-full-v1 --audit results/phylogeny/mafft-audit-full-v1 --output results/phylogeny/mafft-matrix-50-v1
python scripts/compare_alignment_correspondence.py --profile results/phylogeny/profile-alignments-full-v1 --profile-audit results/phylogeny/profile-alignment-audit-full-v1 --mafft results/phylogeny/alignments-full-v1 --mafft-audit results/phylogeny/mafft-audit-full-v1 --output results/phylogeny/alignment-correspondence-v1
python scripts/plot_alignment_correspondence.py --comparison results/phylogeny/alignment-correspondence-v1
```

Correspondence is evaluated on the identical protein residues retained by both 50%-occupancy masks. Each method defines a set of cross-taxon residue pairs occupying the same alignment column. A contingency table of profile-column/MAFFT-column memberships counts shared edges as the sum of n(n−1)/2 over cells, avoiding explicit pair expansion. Edge Jaccard is intersection divided by union; fractions of each method's edges recovered by the other are reported separately. Undefined values (no eligible edges) remain missing. Counts of residues retained by each method and by both distinguish coverage differences from correspondence differences. Per-profile-column summaries support later alignment-sensitive analyses.

Input alignment and Stockholm hashes are checked; complete audit gates, identical source-protein hashes and ungapped residues, profile-to-Stockholm mapping, identical taxon sets and independently recomputed occupancy masks are required. Tests compare the contingency calculation against explicit cross-taxon edge enumeration and check invariance to column renumbering. Agreement is not ground-truth accuracy: both methods can make the same error. This comparison does not itself establish tree robustness or justify a new filtering threshold; topology/model/taxon sensitivities remain required. Plotting uses the existing structural-comparisons environment.

## Marker-tree support audits

`audit_marker_tree_support.py` verifies inference success, matrix/input/tree checksums, exact taxon membership, finite nonnegative branch lengths, and the expected n−3 internal edges of a fully resolved unrooted IQ-TREE representation. Internal bipartitions receive deterministic identities from their canonical taxon side and full taxon universe, independent of the arbitrary trifurcating display root. These identities cannot be compared across different taxon universes without an explicit shared-taxon projection.

Support is SH-aLRT from 1,000 replicates, not bootstrap support or posterior probability. Unreported internal support remains missing; it is never converted to zero. Summaries distinguish total edges, edges with support and edges without support. Counts at 50, 80 and 95 are descriptive thresholds, not a calibrated confidence statement. Rooting, branch dating, weak-edge collapse and species-tree reconciliation are separate pending analyses.

```bash
python scripts/audit_marker_tree_support.py --trees results/phylogeny/marker-gene-trees-v2 --output results/phylogeny/marker-support-snapshot-v1 --allow-incomplete
```

The default requires completion of all expected markers. `--allow-incomplete` produces a diagnostic snapshot with explicit pending marker IDs, never a final species-tree dataset. The first snapshot covers 3/125 markers: 1,422 internal branches, 1,390 with reported support and 32 without. IQ-TREE logs also record identical input sequences; missing-support branches must be assessed in that context without assuming all missing values share one cause. Full branch records remain outside Git, while summary tables and receipts are versioned. Tests verify root-representation invariance, taxon-universe identity, missing-versus-zero support and invalid branch/support rejection.
