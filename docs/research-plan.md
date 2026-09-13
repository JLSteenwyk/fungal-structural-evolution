# Research plan

## Sampling and hypotheses
Target 500 distinct fungal species plus 25 outgroup species. Reserve broad coverage (~300), dense comparisons around replicated transitions (~150), and gap filling (~50); these overlap biologically but the final manifest counts unique species once. Choose genomes using taxonomy, phylogenetic distinctiveness, annotation evidence, assembly quality, and ecological sources. Do not substitute the top 500 quality scores for phylogenetic sampling.

H1: particular branches show excess structural change conditional on sequence divergence and family constraints. H2: duplication permits asymmetric structural divergence. H3: replicated ecological transitions associate with structural change. H4: local functional surfaces change despite conserved global folds. Treat these as hypotheses, not expected findings.

## Analysis sequence and gates
1. Inventory public assemblies, annotations, proteomes and existing structures. Freeze accessions and checksums. Resolve fungal circumscription and close-relative availability; select 25 outgroups with evidence. Assess lineage coverage before finalizing taxa. Missing ecology/quality remains unknown.
2. Obtain selected genome, protein, coding sequence and annotation files. Verify checksums and relationships among IDs. Standardize isoforms without discarding paralogs; check contamination, completeness and gene-model errors. Document quality exceptions for sparse lineages.
3. Infer sequence-based species tree from vetted markers. Compare concatenation and gene-tree approaches; assess support, discordance, compositional heterogeneity, missingness and long-branch sensitivity. Evaluate dating with calibration and topology sensitivity; otherwise report relative changes.
4. Infer protein families and domains, gene trees and duplication/loss reconciliation. Structural clusters do not establish orthology. Use different supported taxon sets for each family.
5. Reuse traceable predictions, predict missing proteins, and record model/database versions, sequence hashes and confidence. Quantify source effects. Analyze domains separately from uncertain inter-domain orientations. Preserve low-confidence proteins in the inventory without forcing them into structural inference.
6. Build protein/domain structural atlas, annotated remote-homology candidates and unknown-family catalog. Benchmark novelty against sequence profiles, structural references, artifacts and sampling gaps.
7. Estimate amino-acid and 3Di branch changes on matched supported genealogies. Benchmark structural-alphabet behavior against direct local/contact comparisons; avoid assuming additive geometric distances. Begin with orthologs, then reconciled duplicates. Account for genealogical uncertainty rather than forcing discordant families onto the species tree.
8. Model branch-specific structural excess with family effects, duration, sampling, length, coverage and prediction uncertainty. Avoid treating pairwise comparisons or shared time denominators as independent evidence. Use null simulations and multiple-testing control; quantify power and false positives within the full dataset.
9. Analyze architecture changes and family turnover separately. Test duplication asymmetry and replicated ecological transitions with phylogenetic control and sensitivity to trait uncertainty. Partition local changes among cores, surfaces, catalytic sites, pockets and supported interfaces.
10. Use codon analyses only where alignment and synonymous divergence support inference. Structural acceleration is not positive selection. Reconstruct ensembles of plausible ancestral amino-acid sequences for selected families; compare predicted structures and propose biochemical validation.

## Required robustness
Sequence-derived structures are not independent of sequences. Leading results must survive alternative predictions/settings, confidence/coverage filters, topology uncertainty and available experimental comparisons. Selection of focal cases after discovery must be labeled exploratory; retain independent contrasts for validation where possible.

## Completion
All eight evolutionary analysis objectives in objective.txt require outputs or explicit scientifically justified applicability determinations backed by data. An atlas alone is not project completion. Final deliverables include provenance, trees, uncertainty estimates, tests, figures/tables, executable methods, mechanistic cases, and experimental proposals. Expensive runs require resource estimates first; new paid charges require authorization.
