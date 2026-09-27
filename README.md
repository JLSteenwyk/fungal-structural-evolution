# Fungal structural evolution

Comparative structural genomics across approximately **500 fungal species plus
25 non-fungal outgroups**. We aim to identify where protein structures change
across the fungal phylogeny, how those changes relate to sequence evolution,
and their associations with duplication, domain architecture and ecology.

## Current checkpoint — 27 September 2026

**The project is not complete.** The scheduled ESMFold prediction batches have
finished; this does not mean every fungal protein has a structure or that the
evolutionary analyses have finished. GPU prediction remains paused. Authorized
CPU analyses and background catalog retrieval continue.

The analyzed cohort contains **501 fungal entries and 25 outgroups** (526
entries). The candidate manifest has 527 entries; *Saccharomyces jurei* is
excluded from analysis because no usable annotated proteome was acquired. Twenty-one uncertain fungal labels and two curated hybrids remain
explicit; these entries are not yet established as 501 unique fungal species.
See the [taxon identity review](docs/taxon-identity-sensitivities.md).

| Component | Verified scope | Remaining work |
|---|---|---|
| Completed local predictions | 25,322 ESMFold models across five cohorts, preserving seven prediction configurations | Remaining marker gaps and broader proteome coverage |
| Combined marker availability | September 27 frozen-cache refresh and full reverse-log check leave 3,639 marker/protein records without cached models, representing 3,626 unique sequences; unchanged from September 25 | Availability precedes confidence filtering and is not proteome-wide coverage |
| Whole-proteome AlphaFold catalog | 1,290,278 models linked to 1,319,513 of 5,815,847 representative proteins (22.69%) | Source-specific coverage excludes local ESMFold; confidence qualification and comparative atlas analyses remain incomplete |
| Complete ESMFold paired inputs | Full input readback passed for 122 markers, 294 taxa and 6,758,598 observed paired AA/3Di cells; all 488 fits and their audit completed | Resampling and evolutionary integration |
| Refreshed AlphaFold marker fits | Verified 125-marker collection: 500 fits from 95 unchanged and 30 refitted markers; 61,893 paired branch rows | Full 75,000-draw paired resampling and subsequent audit running; calibrated evolutionary tests remain |
| Direct geometry and fitted tree paths | 2,527,033 pairs across 122 markers; all 10,108,132 tree-path values and 1,464 descriptive rank rows checked | Geometry numerically sampled one pair per marker; phylogenetic dependence, uncertainty and biological acceleration tests remain |
| Candidate domains | All 1,078,592 intervals extracted and atom-audited; database sequence/coordinate readback and 70,537-cluster partition verified; all 594,797 boundary pairs compared; complete family/taxon source joins audited | Confidence and clustering-parameter sensitivity, phylogenetic integration and evolutionary tests |
| Functional correspondences | All 17,105 annotation rows checked separately for both sources: 6,444 observed AlphaFold rows (314 taxa/23 markers) and 5,576 ESMFold rows (283 taxa/21 markers). Recovered AlphaFold accessibility/annotation join verified across all 47,529 sites | Branch/site tests, matched backgrounds and biological interpretation; annotation rows are not independent events |
| Functional-site predictor sensitivity | 150 exact protein coordinates observed in both sources; all partner contexts and 600 comparisons on identical coordinate subsets checked | Native descriptors, side-chain/pocket evidence and broader controls; agreeing predictions are not experimental validation |
| ESMFold accessibility and site coupling | All 25,322 models merged; paired projection verified for 6,758,598 observations. Conditional coupling completed for 122 markers/44,198 sites, including marker resampling and copy-omission sensitivity | Phylogenetic, prediction, alignment and model uncertainty; accessibility does not establish binding interfaces |
| Recovered AlphaFold accessibility | 30,618 models/15,960,692 residues; 9,453,757 paired observations projected and normalized. Full site-summary readback covers 125 markers/47,529 sites | Rate optimization and its audits precede rate/exposure integration and conditional coupling |
| Duplication comparisons | Matched domain measurements and 82,944 record/family/taxon summary rows verified across all 192 settings; full input-order and weighting sensitivity checked | All 412,800 primary whole-protein tasks processed; strict audit stopped on a two-residue RMSD discrepancy. Full diagnosis confirmed 27 two-residue discrepancies; geometry qualification and background comparisons continue. Phylogenetically adjusted effects, uncertainty and biological interpretation remain |

The latest [matched contrast sensitivity](docs/matched-record-sensitivity-20260927.md)
shows that weighting can reverse the descriptive difference in some settings.
[Taxon and reuse inputs](docs/selected-taxon-inputs-20260927.md) preserve the
dependence information used by the [full comparative-model run](docs/full-matched-working-models-20260927.md).
Its verified inventory contains 28,808 unique record inputs and 144,040 tree
fits. Production is running; numerical output checks, full analytic convergence
assessment, unresolved-fit refinement and inferential calibration remain pending.
Partial fits are not final adjusted effects.

The [ordinary-ML nonlinear comparison](docs/matched-ordinary-likelihood-20260927.md)
is now running across 432,120 fits: linear, quadratic and cubic identity terms
on the same observations under five trees. Full output audit and linked
likelihood comparison are queued. These sensitivity fits do not yet establish
model preference, significance or calibrated uncertainty.

The completed [joint-covariate support check](docs/joint-covariate-support-20260927.md)
now covers all 28,808 unique inputs and 82,944 original settings. Twenty-four
tiny negative-weight cases passed a separately recorded nonnegative construction;
original solver flags remain preserved. Numerical overlap does not establish
model adequacy or calibrated uncertainty. A complete reporting export is queued
after analytic validation, retaining all 414,720 setting/tree dispositions.

[Ecological edge/coverage integration](docs/ecology-optimal-edge-states-20260927.md)
finds no matched predictor marker coverage for the three currently mapped focal
taxa. The separate [wood-decay diagnostic](docs/wood-decay-phylogenetic-diagnostic-20260927.md)
retains eight reviewed taxa and uncertain classes; minimum-change locations
depend on coding. Mapping across all 2,000 bootstrap trees is running, followed
by full independent verification. These diagnostics do not establish replicated
origins or ecological structural effects.

Coverage sources and exact limitations:
[completed prediction inventory](docs/completed-prediction-inventory-20260922.md),
[remaining marker gaps](docs/remaining-marker-model-gaps-20260922.md),
[September 27 residual inventory](metadata/current_marker_residual_gap_summary_20260927.json),
[September 26 completion evidence](docs/recovery-20260926.md),
[whole-proteome catalog](docs/whole-proteome-structure-catalog-20260922.md),
[domain registry and clustering](docs/whole-proteome-structure-domain-registry-20260923.md),
[functional sites](docs/functional-site-workflow.md), and
[ESMFold accessibility](docs/residue-accessibility.md), and
[recovered AlphaFold accessibility](docs/recovered-afdb-accessibility-20260927.md).

Supported species-tree comparisons and branch resampling remain active.
Completed ESMFold conditional site-rate coupling supports a positive sequence-rate
association, including after omission of one copy-flagged marker. This is a
conditional association, with model-dependent interaction evidence, not a causal
effect or a physical displacement estimate. See the
[full-cohort coupling results](docs/conditional-site-coupling.md).
The full path and rank benchmark is descriptive, with shared ancestry and
prediction uncertainty still requiring treatment. See the
[September 26 checkpoint](docs/recovery-20260926.md) for completed stages and
[orthogroup validation](docs/hierarchical-orthogroup-validation-20260922.md)
for reconciliation limitations. Queued stages require successful, source-bound
completion records before advancing.

## Evolutionary objectives still open

| Objective | Work still required before completion |
|---|---|
| Branches and clades with elevated structural change | Complete branch uncertainty, direct-coordinate benchmarking, calibrated tests and multiple-testing control |
| Sequence–structure coupling | Extend full ESMFold conditional results to the AlphaFold cohort; propagate topology, alignment and prediction uncertainty |
| Domain evolution | Complete structural comparisons and reconcile gains, losses, fusions and rearrangements separately from family turnover |
| Duplication-associated divergence | Complete protein/domain comparisons, model/reference quality sensitivity, matched controls and phylogenetically controlled divergence/asymmetry tests |
| Ecological transitions | Establish usable replicated contrasts, propagate trait uncertainty and test phylogenetically controlled associations |
| Functional locations | Integrate surface/core and functional-site inputs with branch/site changes; assess supported pockets and interfaces where evidence permits |
| Selection | Resolve codon eligibility, saturation, copy and optimization concerns before appropriate tests; the current 1,632 local codon cases have all-case optimization checks running and do not yet establish eligibility |
| Ancestral and mechanistic cases | Select supported exploratory cases, propagate ancestral uncertainty, obtain authorized predictions and formulate testable functional hypotheses |

The [original objective](docs/objective.txt) and [research plan](docs/research-plan.md)
retain the full scope. Earlier subset analyses provide intermediate evidence;
they do not substitute for these requirements. Detailed chronology and
completion receipts are indexed in [progress](docs/progress.md).

## Project records

- [Original objective](docs/objective.txt)
- [Research plan](docs/research-plan.md)
- [Progress and completion evidence](docs/progress.md)
- [Open scientific milestones](docs/scientific-milestones-20260926.md)
- [Duplication coverage and ascertainment](docs/duplication-sampling-coverage-20260926.md)
- [Duplicate/reference domain controls](docs/duplication-domain-controls-20260926.md)
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

## Figures and earlier cohort-specific checkpoints

The figures and analysis summaries below retain their original dataset scopes.
Some describe earlier execution states; use the current checkpoint above and
linked workflow completion records for current status.

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
