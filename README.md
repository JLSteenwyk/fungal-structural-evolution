# Fungal structural evolution

Comparative structural genomics of approximately **500 fungal species plus 25 non-fungal outgroups**. The central question is where structural evolution accelerates or decouples from sequence evolution, and how these changes relate to duplication and ecology.

Status: proteomes and broad-marker QC are complete for **501 fungal entries and 25 genome-backed outgroups (526 entries)**. A [taxon-label audit](metadata/taxon_label_review.tsv) flags one explicit hybrid and 21 incompletely identified fungal entries for species-level review. All 125 profile and MAFFT marker alignments and their coding-sequence projections are complete. Both full-taxon guide-tree searches, supported marker trees, full orthology assignment and whole-proteome domain searches are running. Supported species trees, reconciled orthology and the full evolutionary analyses remain incomplete. Full-scale design; no separate pilot.

Confidence-qualified structural coverage now reaches **458 entries (440 fungal entries and 18 outgroups)** across the separate AlphaFold and local ESMFold snapshots. The local 5,121-model snapshot adds 4,669 usable taxon–marker combinations and 136 previously unrepresented entries; prediction queues continue on both GPUs. These are coverage counts, not fully sampled or supported evolutionary results. [Coverage evidence and source comparison](docs/structure-prediction-workflow.md#local-confidence-qualification-and-coverage-contribution).

A completed 266-protein predictor control found **14.3% median structural-state disagreement** among 187 confidence-qualified proteins with identical amino-acid sequences. This source effect must be assessed before interpreting branch differences. Supported sequence-tree and structural-model fits are running separately for 124 AlphaFold and 72 local marker alignments. [Predictor control methods and limitations](docs/structure-prediction-workflow.md#same-sequence-structural-alphabet-predictor-control).

Assembly-matched CDS acquisition and strict translation auditing are complete for all 519 NCBI-backed taxa (5,752,457 exact translations). Separate genome/translation audits cover the seven external taxa, retaining exceptions and partial coding boundaries; see the [coding-sequence workflow](docs/coding-sequence-workflow.md). All 125 marker codon alignments have passed independent translation readback; selection analyses remain pending.

## Project records
- [Original objective](docs/objective.txt)
- [Research plan](docs/research-plan.md)
- [Progress and completion evidence](docs/progress.md)
- [Decisions and questions](docs/decisions.md)
- [Literature](docs/bibliography.md)
- [Resource assessment](docs/resources.md)
- [Phylogenetic workflow](docs/phylogenetic-workflow.md)
- [Gene and isoform reconciliation](docs/gene-isoform-mapping.md)
- [Assembly-quality workflow](docs/assembly-quality-workflow.md)
- [Orthology workflow](docs/orthology-workflow.md)
- [Domain annotation workflow](docs/domain-annotation-workflow.md)
- [Domain geometry and placement confidence](docs/domain-structure-comparisons.md)
- [Structural alphabet benchmark](docs/structural-alphabet-benchmark.md)
- [Fitted tree paths and direct geometry](docs/tree-path-geometry.md)
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

The expanded source-specific catalog now contains 13,153 GDM models linked to 13,654 marker proteins in 322 taxa. A [full-design coverage audit](docs/prediction-source-controls.md#expanded-model-catalog-and-full-design-coverage) retains all 526 taxa and separates acquisition candidates, query gaps and missing marker sequences. Expanded residue mapping is complete with 4,895,692 residue links; native structural-alphabet extraction, export checks and coordinate-feature reconstruction are complete; PAE retrieval, final mapping-bound validation and per-model confidence qualification are complete for all 13,153 models, retaining 4,714,151 native states under joint six-residue pLDDT ≥70 and directional PAE ≤10 filters. Expanded paired alignments are complete for 124 markers (5–135 taxa each), retaining 304 fungal entries and 18 outgroups across 23 lineage/role groups. The supported sequence-topology/structural-branch fitting batch is running. Four groups still have no usable paired coverage; evolutionary results remain pending.

A [full-sampling gene-copy annotation audit](docs/gene-isoform-mapping.md#full-sampling-busco-copy-annotation-audit) separates 312 raw duplicated BUSCO calls attributable to one annotated gene’s multiple products from 1,989 calls spanning multiple annotated genes and four unresolved calls. Biological duplication, assembly redundancy and contamination remain to be assessed.

Both full-taxon alignment strategies now have guide-tree searches running. A [taxon-coverage sensitivity audit](docs/phylogenetic-workflow.md#full-sampling-taxon-coverage-sensitivity-definitions) shows that a 50% coverage filter would eliminate all sampled Microsporidia and Olpidiomycota; the production design retains all 526 taxa. Supported trees and filtered-tree comparisons remain pending.

[Assembly-matched coding-sequence acquisition and translation auditing](docs/coding-sequence-workflow.md) are complete for the 519 NCBI-backed taxa, with separately audited external CDS sources and explicit exceptions. All 125 marker codon alignments passed translation readback. Selection modeling still requires orthology, alignment and divergence review; the completed gene-tree subset already flags conflicts with some genus-based grouping assumptions.
