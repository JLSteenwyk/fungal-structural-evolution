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

## Full-matrix alignment sensitivity tree execution

The full **526-taxon, 63,750-column MAFFT matrix** has now entered its own IQ-TREE guide search, using the same LG+F+G4 model, seed 20260913, 16 threads and 32 GB memory limit as the profile-based guide. All artifacts in both matrix receipts passed checksum validation; their exact 526 unique taxa agree and every sequence has its declared matrix length. This independent tree search advances alignment-method sensitivity across the complete design.

```bash
python scripts/run_initial_species_tree.py \
  --matrix results/phylogeny/mafft-matrix-50-v1 \
  --output results/phylogeny/initial-mafft-guide-v1
```

The prelaunch estimate allows 4–24 hours and 5 GB output space on the existing host, with uncertain tree-search runtime. IQ-TREE's actual initialization estimates 20,661 MB RAM, within the limit, and has entered likelihood optimization. It reports 31 sequences with more than 50% gaps/ambiguity and nominal composition-test failures for all 526 sequences. Those flags reinforce the requirement for compositional/model sensitivity; they are not automatic taxon exclusions or calibrated evidence that a particular topology is wrong. The homogeneous-model guide has no bootstrap support and is not a final species tree. Configurations and resource evidence: `metadata/mafft_guide_run_config.json` and `metadata/mafft_guide_resource_plan.json`. The guide's terminal receipt and subsequent topology/support comparison remain pending.

## Full-sampling taxon-coverage sensitivity definitions

A direct readback of both matrices reproduces the recorded unambiguous amino-acid counts for every taxon. The audit joins the raw single-copy, duplicated, fragmented and missing marker counts so that sparse matrix occupancy is not equated with poor overall genome completeness. The threshold is applied to the fraction of unambiguous sites in **both** matrices; these matrices have different alignment scopes and denominators.

| Minimum coverage in both matrices | Fungal entries retained | Outgroups retained | Manifest groups entirely lost |
| --- | ---: | ---: | --- |
| 10% | 499 | 23 | None |
| 30% | 484 | 23 | Olpidiomycota |
| 50% | 472 | 23 | Microsporidia; Olpidiomycota |
| 70% | 452 | 20 | Microsporidia; Olpidiomycota; Sanchytriomycota |

These are explicit sensitivity definitions, **not calibrated reliability thresholds or final filtering decisions**. The production matrices and ongoing full-taxon guides retain all 526 taxa. The table shows why more aggressive filtering cannot silently replace the broad sampling design. Manifest bins have different taxonomic ranks.

Four taxa fall below 10% in both alignments: Amoeboaphelidium protococcarum F1243177 (164 profile / 204 MAFFT unambiguous sites; one single-copy and 119 duplicated markers), Pirum gemmata OFS5426506 (1,224 / 1,715 sites; seven single-copy markers), Abeoforma whisleri OFS5426458 (1,804 / 2,272; seven single-copy markers), and Cryoendolithus antarcticus F1507870 (4,203 / 5,765; 15 single-copy and 110 duplicated markers). Presence in a tree does not establish reliable placement. The high-duplication taxa need gene-copy/assembly review; the sparse outgroups need marker/placement sensitivity and appropriate family-specific inclusion.

```bash
python scripts/assess_taxon_phylogenetic_coverage.py \
  --profile results/phylogeny/profile-matrix-50-v1 \
  --mafft results/phylogeny/mafft-matrix-50-v1 \
  --copies results/qc/busco-gene-copies-v1 \
  --output results/phylogeny/taxon-coverage-sensitivity-v1
```

All source artifacts, exact taxon membership, matrix length/alphabet, per-taxon counts and output hashes were checked. The output includes a complete 526-taxon table and explicit threshold membership for each taxon. No filtered tree has yet been inferred, so these tables alone do not establish topology sensitivity. Versioned summaries are `metadata/taxon_matrix_coverage.tsv`, `metadata/taxon_filter_lineage_sensitivity.tsv` and `metadata/taxon_matrix_coverage_sensitivity_receipt.json`.

## Expanded marker-tree support checkpoint

The second immutable support snapshot validates 12 completed marker trees and 5,742 internal branches. There are 5,630 branches with reported SH-aLRT support and 112 with no reported value; 3,534 reach SH-aLRT 80. Missing support remains distinct from zero, and SH-aLRT is not a bootstrap percentage. All source tree/alignment identities and branch invariants passed the audit; output artifacts passed independent checksum readback.

```bash
python scripts/audit_marker_tree_support.py \
  --trees results/phylogeny/marker-gene-trees-v2 \
  --output results/phylogeny/marker-support-snapshot-v2 --allow-incomplete
```

The 113 pending marker IDs remain explicit in the receipt. Trees that finish early need not be representative of the full marker set. This checkpoint does not establish species-tree support, gene concordance factors, reconciliation or the cause of any discordance. The earlier three-tree snapshot remains preserved. New versioned records: `metadata/marker_tree_support_snapshot_v2.tsv` and `metadata/marker_tree_support_snapshot_v2_receipt.json`.


## Full-matrix mixture sensitivity: resource assessment

Both 526-taxon guide searches remain active. Before the next phylogenetic stage,
`scripts/prepare_species_mixture_resources.py` verified both complete matrix
artifact manifests, equal tip sets and alignment lengths, and froze their guide
startup resource reports. The profile matrix has 48,817 distinct patterns and a
reported homogeneous-model estimate of 15,875 MB; MAFFT has 63,543 patterns and
20,661 MB. These are program estimates, not measured resident memory.

A linear component-count scenario, including the empirical profile, gives:

| Matrix | C20 scenario MB | C20 with 25% allowance MB | C60 scenario MB |
|---|---:|---:|---:|
| Profile | 333,375 | 416,718.75 | 968,375 |
| MAFFT | 433,881 | 542,351.25 | 1,260,321 |

Scaling is a planning assumption, not a guaranteed upper bound. The next feasible
sensitivity is C20-PMSF on both full matrices, crossed with both completed guide
topologies to expose dependence on guide choice. IQ-TREE derives site profiles
from a mixture model and guide tree; its first PMSF phase retains mixture memory
requirements. [IQ-TREE documentation](https://iqtree.github.io/doc/Complex-Models#site-specific-frequency-models);
[Wang et al. 2018](https://doi.org/10.1093/sysbio/syx068).

Plan one run at a time, 16 threads, a 600G IQ-TREE memory limit and at least
750 GiB available host memory before launch, with 100 GB output allowance.
Treating the reported MB conservatively as MiB, both C20 scenarios plus allowance
fit this limit. The initial v1 plan used 512G; arithmetic review showed that the
MAFFT allowance could exceed it, so v2 supersedes it before any launch. A broad
24–336 hour planning range per supported run is unmeasured. Unrestricted C60
requires revisiting memory strategy; reducing taxa is not the chosen remedy.

The resource preparation and independent arithmetic/hash readback passed.
Receipt: `metadata/species_mixture_resource_receipt.json`; frozen logs:
`results/phylogeny/species-mixture-resources-v2`. No mixture fit has launched.
Completed guide receipts and frozen, validated trees are prerequisites. This
sensitivity does not replace partitioned analysis, gene concordance, across-
lineage composition assessment or taxon/marker sensitivity.


## Profile guide completed and first full-matrix PMSF run launched

The 49,027-site profile guide completed after 55,987.5 seconds of wrapper time.
`scripts/audit_species_guide.py` verified the full input manifest, all 526 tips
(501 fungal entries and 25 outgroups), 1,049 finite nonnegative edges, exact
report/tree edge agreement and total branch length 146.2086028457. The reported
log likelihood is −16,717,819.5477 under LG+F+G4. Nominal composition failures
are 520/526; the frozen taxon table preserves the six passing entries as well.
It is unrooted and has no support estimates. Receipts:
`metadata/profile_guide_execution_receipt.json` and
`metadata/profile_guide_full_audit_receipt.json`.

The completed guide now supplies the first of four planned full-matrix C20-PMSF
sensitivity combinations: profile alignment with profile guide. The earlier
both-guides prerequisite is needed for completing the crossed comparison;
this first combination can begin with its own completed, audited guide.
`scripts/run_species_pmsf.py` validated input hashes, the exact taxon grid,
current memory and disk headroom and the pinned IQ-TREE 3.0.1 executable.
It requests 16 threads, a 600G memory limit, 1,000 SH-aLRT and 1,000 UFB
replicates with bootstrap NNI and saved bootstrap trees. A shared execution lock
prevents concurrent full-matrix PMSF runs through this script. IQ-TREE entered
site-profile estimation and reported 333,579 MB required, consistent with the
resource scenario; this is a program estimate, not measured peak memory.

Launch config: `metadata/pmsf_profile_profile_launch_config.json`.
Live output: `results/phylogeny/pmsf-profile-profile-v1`.
PMSF support conditions on estimated site profiles. Full output/model/support
validation remains necessary, as do the other guide/alignment combinations,
partitioned analyses, discordance and across-lineage composition sensitivities.
No supported final species tree is claimed.

```bash
python scripts/audit_species_guide.py \
  --guide results/phylogeny/initial-guide-v1 \
  --matrix results/phylogeny/profile-matrix-50-v1 \
  --output results/phylogeny/profile-guide-audit-v1
python scripts/run_species_pmsf.py \
  --matrix results/phylogeny/profile-matrix-50-v1 \
  --guide-audit results/phylogeny/profile-guide-audit-v1 \
  --resources results/phylogeny/species-mixture-resources-v2 \
  --output results/phylogeny/pmsf-profile-profile-v1
```

## Second full-taxon guide completed; crossed PMSF continuation queued

The MAFFT-based guide completed in 57,023.9 seconds. Its full readback audit
verifies the 526-taxon manifest (501 fungal entries plus 25 outgroups), 63,750
alignment columns, 1,049 finite nonnegative branches and exact report/tree
agreement. Total branch length is 145.6696520524 under LG+F+G4; reported log
likelihood is −22,762,350.7173. All 526 sequences fail the nominal composition
screen. This is an unsupported unrooted guide, not a final species phylogeny.
Its likelihood must not be compared directly with the different profile alignment.

The two guides share 470 of their 523 internal splits. Each contains 53 splits
absent from the other: unrooted Robinson–Foulds distance 106, normalized by their
1,046 total internal splits to 0.10134. Independent graph traversal verifies all
2,098 edge bipartitions. Across 138,075 taxon pairs, tree-path rank correlation is
0.98911; 200 deterministic path calculations were independently checked with
Biopython traversal. This is a descriptive comparison, not independent-pair
inference, support or a test identifying which guide is correct. Differences can
reflect alignment, search and model effects.

For the FCS-flagged N. cerealis entry, the three nearest taxa by tree-path length
are other Naganishia entries in both guides. This concatenated-guide observation
does not resolve the provenance of individual flagged proteins or rule out a
mixed-source assembly. The 29 affected codon cases and 16 earlier ESMFold marker
observations still require their documented source-quality review.

```bash
python scripts/audit_species_guide.py --guide results/phylogeny/initial-mafft-guide-v1 --matrix results/phylogeny/mafft-matrix-50-v1 --output results/phylogeny/mafft-guide-audit-v1
python scripts/compare_full_species_guides.py --profile results/phylogeny/profile-guide-audit-v1 --mafft results/phylogeny/mafft-guide-audit-v1 --manifest metadata/analysis_manifest.tsv --output results/phylogeny/full-guide-comparison-v1
```

Both audited guides now satisfy the prerequisite for the previously planned
2×2 alignment/guide PMSF comparison. A verified live controller waits for the
original profile-alignment/profile-guide producer and requires its successful
completion receipt. It then runs profile/MAFFT, MAFFT/profile and MAFFT/MAFFT
serially. Each uses the existing runner's exclusive lock, 16 threads, 600G limit,
750 GiB minimum available-memory check, 100 GB output allowance and 24–336 hour
planning window per run. These are conservative planning scenarios, not runtime
promises. No concurrent full-matrix PMSF job is added.

```bash
python scripts/advance_crossed_species_pmsf.py --config metadata/crossed_species_pmsf_controller_config.json
```

The controller is already running; do not launch another copy. It pins both
matrix/guide artifact sets, scripts and the first run configuration, validates
successful stage receipts and checkpoints each completed stage. A failed or
uncheckpointed partial stage requires review. Full model/profile/support audits
remain necessary after execution. This is the baseline full-manifest crossed
comparison; separate FCS and taxon-identity sensitivities remain required.

### Completed-marker support snapshot v4

The fourth immutable support snapshot audits 48 of the 125 planned individual
marker trees. It contains 22,705 internal splits: 22,454 with reported SH-aLRT
support and 251 with no reported support value. Missing support is retained as
missing and must not be interpreted as zero. These are 1,000-replicate SH-aLRT
values, not bootstrap percentages or posterior probabilities.

An additional independent readback reconstructed every retained sequence from
the original profile matrix under the established coverage rule, checked all
excluded-taxon counts, and verified 8,465,734 retained alignment characters.
It recovered every internal bipartition by deleting an edge from an undirected
tree graph and traversing the resulting component, then checked split identity,
attached support and branch length against the snapshot. Fitted model labels
were checked against the allowed LG/WAG/JTT+F+G4 family; report and coverage-file
hashes are retained for every completed marker.

```bash
python scripts/audit_marker_tree_support.py --trees results/phylogeny/marker-gene-trees-v2 --output results/phylogeny/marker-tree-support-snapshot-v4 --allow-incomplete
python scripts/readback_marker_tree_snapshot.py --snapshot results/phylogeny/marker-tree-support-snapshot-v4 --trees results/phylogeny/marker-gene-trees-v2 --matrix results/phylogeny/profile-matrix-50-v1 --output results/phylogeny/marker-tree-support-readback-v4
```

Receipts and marker summaries are versioned in
`metadata/marker_tree_support_{snapshot,readback}_v4*`; full split tables remain
outside Git. The 77 pending marker trees are not failed or excluded trees.
Completion order can depend on alignment length, taxon coverage and inference
difficulty, so this subset cannot stand in for full-batch discordance estimates.
The checks establish input and split integrity, not model adequacy, rooting,
orthology, species-tree concordance or reconciliation. FCS-flagged observations
remain present in these original baseline trees and require their planned
sensitivity and copy reviews.
