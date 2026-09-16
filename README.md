# Fungal structural evolution

Comparative structural genomics of approximately **500 fungal species plus 25 non-fungal outgroups**. The central question is where structural evolution accelerates or decouples from sequence evolution, and how these changes relate to duplication and ecology.

Status: proteomes and broad-marker QC are complete for **501 fungal entries and 25 genome-backed outgroups (526 entries)**. The [taxon-label audit](metadata/taxon_label_review.tsv) flags 21 incompletely identified fungal entries; [assembly-linked literature review](docs/taxon-identity-sensitivities.md) identifies two curated hybrids. These entries are not yet established as 501 unique fungal species. All 125 profile and MAFFT marker alignments and their coding-sequence projections are complete. Both full-taxon profile and MAFFT guide trees are complete and audited. The first supported C20-PMSF species-tree analysis is running, with three crossed alignment/guide analyses queued. Supported marker trees, full orthology assignment and whole-proteome domain searches remain running; the latest audited marker-tree snapshot covers 70/125 markers. Supported species trees, reconciled orthology and the full evolutionary analyses remain incomplete. Full-scale design; no separate pilot.

Confidence-qualified paired-marker coverage reaches **525 entries (500 fungal entries and all 25 outgroups)** across separate AlphaFold and ESMFold datasets, with 22,677 usable taxon–marker combinations in their descriptive union. The three audited ESMFold acquisition cohorts are now integrated: 10,048 models yield 89 ready markers across 292 taxa. Paired branch fits and direct coordinate comparisons are audited; expanded ESMFold site-rate optimization, 89-marker conditional coupling and an 88-marker copy-review omission sensitivity are complete, including marker resampling. AlphaFold optimization and its reviewed coupling handoffs remain in progress. *Amoeboaphelidium protococcarum* remains uncovered. These counts do not establish 500 unique species, uniform coverage or pooled evolutionary inference. Prediction-source controls remain necessary. [Coverage records](metadata/combined_predictor_coverage_union.json).

A completed 266-protein predictor control found **14.3% median structural-state disagreement** among 187 confidence-qualified proteins with identical amino-acid sequences. This source effect must be assessed before interpreting branch differences. Fits for the earlier 72-marker local cohort are complete, with conditional resampling and direct-geometry comparisons; the separate 124-marker AlphaFold batch has completed all 496 fits, its full output audit and direct-geometry benchmark. Paired site/block resampling is now running for all 124 AlphaFold markers. The integrated 89-marker ESMFold cohort is also undergoing supported branch inference. [Predictor control methods and limitations](docs/structure-prediction-workflow.md#same-sequence-structural-alphabet-predictor-control).

The expanded ESMFold analysis now includes 24 models across 22,205 sites and a matched 88-marker/21,904-site omission sensitivity. Both retain positive sequence-rate coefficients across all specifications, with unadjusted marker-bootstrap intervals above zero; interaction intervals include zero. These remain conditional associations, and copy reconciliation is unresolved. The earlier 72-marker snapshot remains available: 24 models across 16,571 sites, with marker intercepts and coverage, confidence and composition controls. Sequence-rate coefficients are positive across specifications, but the full models explain only 5–7% of within-marker variation. These are exploratory associations conditional on fixed trees and estimated rates. Omitting FCS-overlapping observations from 16 markers preserves positive sequence-rate coefficients in all 24 models, but 20 meet the existing BH threshold versus 23 at baseline. This does not test coefficient differences; uncertainty propagation and prediction circularity remain unresolved. [Methods, results and figure](docs/conditional-site-coupling.md).

[Assembly contamination report review](docs/assembly-quality-workflow.md) now includes 518 usable NCBI FCS reports, one unresolved publisher-checksum mismatch and seven external taxa without these NCBI reports. Exact CDS mapping identifies 53 marker observations overlapping EXCLUDE/FIX/TRIM regions in two fungal assemblies, plus one REVIEW-only observation. This does not by itself confirm contamination. Paired-marker omission refits are complete for both ESMFold cohorts; alternative-copy review trees and full-taxon guide sensitivities are running. Retained-site geometry, accessibility and exposure-change checks are complete for the affected paired markers; revised rate and coupling fits are complete for the earlier 72-marker cohort, while expanded baseline coupling and copy-review omission are complete; expanded FCS-specific coupling remains outstanding.

Assembly-matched CDS acquisition and strict translation auditing are complete for all 519 NCBI-backed taxa (5,752,457 exact translations). Separate genome/translation audits cover the seven external taxa, retaining exceptions and partial coding boundaries; see the [coding-sequence workflow](docs/coding-sequence-workflow.md). All 125 marker codon alignments have passed independent translation readback. Within-genus diagnostics now include 1,655 supported nucleotide trees and 1,655 global MG94 fits, with full bootstrap and saved-likelihood readback plus independently checked site-normalized branch distances. All 13,240 nuisance-reoptimized branch-parameter profile fits and saved-model readbacks are complete and audited; four cases require further optimization review. Twenty-nine codon cases include FCS-flagged observations and fall below the existing four-taxon minimum if those observations are omitted. Saturation, gene-copy and other eligibility reviews remain pending before selection tests.

The experimental-reference benchmark now includes **79 distinct predicted proteins**, with 71/63/18 eligible at joint confidence cutoffs 0/70/90. The [updated benchmark figure](docs/figures/combined_experimental_controls.svg) shows all four length tiers, exclusions and the same-protein sensitivity; [methods and limitations](docs/experimental-structure-workflow.md) explain why these selected references do not establish unbiased predictive accuracy.

## Project records

- [Original objective](docs/objective.txt)
- [Research plan](docs/research-plan.md)
- [Progress and completion evidence](docs/progress.md)
- [Decisions and questions](docs/decisions.md)
- [Literature](docs/bibliography.md)
- [Resource assessment](docs/resources.md)
- [Phylogenetic workflow](docs/phylogenetic-workflow.md)
- [Dating evidence and calibration review](docs/dating-workflow.md)
- [Gene and isoform reconciliation](docs/gene-isoform-mapping.md)
- [Assembly-quality workflow](docs/assembly-quality-workflow.md)
- [Orthology workflow](docs/orthology-workflow.md)
- [Domain annotation workflow](docs/domain-annotation-workflow.md)
- [Domain geometry and placement confidence](docs/domain-structure-comparisons.md)
- [Structural alphabet benchmark](docs/structural-alphabet-benchmark.md)
- [Fitted tree paths and direct geometry](docs/tree-path-geometry.md)
- [Conditional site-level coupling](docs/conditional-site-coupling.md)
- [Prediction-source controls](docs/prediction-source-controls.md)
- [Final outgroup query resolution](docs/chromosphaera-taxonomy-resolution.md)
- [Local structure prediction](docs/structure-prediction-workflow.md)
- [Ecological evidence and curation](docs/ecology-evidence-workflow.md)
- [Structure-to-phylogeny residue mapping](docs/marker-structure-integration.md)
- [Methods draft: executed work and pending analyses](docs/methods-draft.md)
- [Metadata definitions](metadata/README.md)

## Reproduction
Python 3.10+ standard library suffices for NCBI metadata inventory; openpyxl is used to import published supplementary tables. Run `make inventory` to retrieve public NCBI fungal catalogs, record checksums, and generate candidate tables. These are discovery catalogs, not a final sample. Raw downloads remain in ignored `data/`; versioned metadata records source URLs and hashes.

BUSCO uses the isolated environment specified in `environments/busco.yml`; marker workflows additionally use Biopython and MAFFT. Commands and validation gates are documented in the phylogenetic workflow. Large files must not enter Git history. GitHub: public repository https://github.com/JLSteenwyk/fungal-structural-evolution (user-confirmed).

## Full-dataset QC

![Sampling and broad marker recovery](docs/figures/marker_recovery.svg)

These are raw-proteome results for the 125-marker broad eukaryotic panel, not a universal assembly-quality score. Reproduce with `python scripts/summarize_busco.py` and `python scripts/plot_full_dataset_qc.py`; the latter writes PNG, SVG, PDF and table artifacts under `results/qc/`. Copy its SVG to `docs/figures/marker_recovery.svg` to refresh the displayed figure. Plot dependencies are in `environments/qc-figures.yml`.

## Exploratory structure confidence assessment

![PAE confidence sensitivity](docs/figures/pae_sensitivity.svg)

Version-matched PAE matrices were validated for 422 models, enabling confidence sensitivity for all 858 marker taxon pairs qualifying at pLDDT ≥70. The figure compares residue-distance changes before and after a directional PAE filter. Filters change which residue pairs are compared; these descriptive measurements do not establish branch rates, domain movement or adaptation. Reproduction and interpretation are in the [structure integration workflow](docs/marker-structure-integration.md).

## Alignment-method sensitivity

![Profile and MAFFT correspondence](docs/figures/alignment_correspondence.svg)

All 125 markers have completed profile and MAFFT alignments. Their alternative 526-taxon matrices retain 49,027 and 63,750 columns, respectively. Median residue-pair Jaccard agreement among residues retained by both methods is 0.945, with substantial disagreement in a few markers. Agreement is not accuracy; tree and structural-result sensitivities remain pending. See the [phylogenetic workflow](docs/phylogenetic-workflow.md).

Domain annotation is now linked to structural geometry for the existing AFDB marker snapshot: 738 comparisons across 84 Pfam domains, plus 375 domain-pair/taxon-pair placement assessments. These measurements separate internal domain geometry from global-fit effects and uncertain interdomain placement; they are not yet evolutionary rate estimates. [Methods, figure and an artifact-control example](docs/domain-structure-comparisons.md).

Coordinate-derived 3Di encodings and all ten native features have been audited for the 422-model comparison snapshot. The matched-site sequence/3Di/geometry benchmark covers the 858 whole-marker and 738 domain comparisons under three confidence regimes. Invalid terminal states are explicitly masked and spatial-partner confidence is retained. [Benchmark figure, methods and limitations](docs/structural-alphabet-benchmark.md); matched-topology point estimates now cover 52 markers and 416 branches under three structural models. [Paired phylogenetic methods and audit](docs/paired-structural-phylogenetics.md); 31,200 paired site/block resampling draws have now been audited, providing conditional interval sensitivities. Nonlocal feature dependence, broader uncertainty, model adequacy and acceleration tests remain unresolved.

The expanded source-specific catalog now contains 13,153 GDM models linked to 13,654 marker proteins in 322 taxa. A [full-design coverage audit](docs/prediction-source-controls.md#expanded-model-catalog-and-full-design-coverage) retains all 526 taxa and separates acquisition candidates, query gaps and missing marker sequences. Expanded residue mapping is complete with 4,895,692 residue links; native structural-alphabet extraction, export checks and coordinate-feature reconstruction are complete; PAE retrieval, final mapping-bound validation and per-model confidence qualification are complete for all 13,153 models, retaining 4,714,151 native states under joint six-residue pLDDT ≥70 and directional PAE ≤10 filters. Expanded paired alignments are complete for 124 markers (5–135 taxa each), retaining 304 fungal entries and 18 outgroups across 23 lineage/role groups. The paired sequence-topology/structural-branch fitting batch is audited; rate-model sensitivity and downstream coupling remain in progress. Four groups still have no usable paired coverage; evolutionary results remain pending.

A [full-sampling gene-copy annotation audit](docs/gene-isoform-mapping.md#full-sampling-busco-copy-annotation-audit) separates 312 raw duplicated BUSCO calls attributable to one annotated gene’s multiple products from 1,989 calls spanning multiple annotated genes and four unresolved calls. Biological duplication, assembly redundancy and contamination remain to be assessed.

Both full-taxon homogeneous guides are complete and audited, sharing 470 of their 523 internal splits. The first supported C20-PMSF analysis and separate full-taxon FCS mask guide sensitivities are running. A [taxon-coverage sensitivity audit](docs/phylogenetic-workflow.md#full-sampling-taxon-coverage-sensitivity-definitions) shows that a 50% coverage filter would eliminate all sampled Microsporidia and Olpidiomycota; the production design retains all 526 taxa. Supported trees and filtered-tree comparisons remain pending.

[Assembly-matched coding-sequence acquisition and translation auditing](docs/coding-sequence-workflow.md) are complete for the 519 NCBI-backed taxa, with separately audited external CDS sources and explicit exceptions. All 125 marker codon alignments passed translation readback. Selection modeling still requires orthology, alignment and divergence review; the completed gene-tree subset already flags conflicts with some genus-based grouping assumptions.

The full local ESMFold resampling run now has an audited 43,200 paired draws
across 72 markers (16 unestimable; 86,368 validated fits). Joint path summaries
are complete for 259,780 accepted structural comparisons and 6,813 exclusions;
independent path readback passed, including 86,368 tree traversals. [Methods and tracked receipts](docs/tree-path-geometry.md).
These conditional sensitivities are an intermediate result. Species phylogeny,
full structural coverage and the evolutionary hypothesis tests remain incomplete.


[Lineage-specific completeness QC](docs/lineage-completeness-workflow.md) is
complete for all 501 fungal entries. The joined quality table retains the 25
outgroups with broad-panel scores; panel differences and missing markers require
lineage-aware interpretation. The [FCS report/CDS overlap audit](docs/assembly-quality-workflow.md) is complete; biological source review and contamination sensitivity analyses remain pending.
