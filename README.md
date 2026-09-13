# Fungal structural evolution

Comparative structural genomics of approximately **500 fungal species plus 25 non-fungal outgroups**. The central question is where structural evolution accelerates or decouples from sequence evolution, and how these changes relate to duplication and ecology.

Status: proteomes acquired and broad-marker QC completed for **501 fungal entries and 25 genome-backed outgroups (526 taxa)**. A [taxon-label audit](metadata/taxon_label_review.tsv) flags one hybrid and 21 incompletely identified fungal entries for species-level review. Gene-representative preparation is complete with unresolved mappings flagged. All 125 profile-based marker alignments are complete; an initial tree search is running on a 49,027-position matrix. OrthoFinder reference-core inference has started with all remaining taxa prepared for assignment. MAFFT sensitivity alignments and existing-structure acquisition continue. Supported species trees, full orthology and evolutionary analyses remain pending. Full-scale design; no separate pilot.

## Project records
- [Original objective](docs/objective.txt)
- [Research plan](docs/research-plan.md)
- [Progress and completion evidence](docs/progress.md)
- [Decisions and questions](docs/decisions.md)
- [Literature](docs/bibliography.md)
- [Resource assessment](docs/resources.md)
- [Phylogenetic workflow](docs/phylogenetic-workflow.md)
- [Gene and isoform reconciliation](docs/gene-isoform-mapping.md)
- [Orthology workflow](docs/orthology-workflow.md)
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
