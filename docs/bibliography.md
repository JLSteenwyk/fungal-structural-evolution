# Annotated bibliography

- Wu et al. 2026. Structural genomics across insects. https://doi.org/10.1038/s41422-026-01220-0 — Atlas plus phylogenetic context and functional validation; pairwise remote-homology counts are not unique proteins.
- Lemke et al. 2025. The role of metabolism in shaping enzyme structures over 400 million years. https://doi.org/10.1038/s41586-025-09205-6 — Yeast enzyme precedent linking structural context with metabolic properties; informs local analyses.
- Derbyshire & Raffaele 2023. https://doi.org/10.1038/s41467-023-40949-9 — Fungal orphan effectors, ancestral reconstruction and surface frustration; relevant to case-study design.
- Seong & Krasileva 2023. https://doi.org/10.1038/s41564-022-01287-6 — Comparative fungal effector structures and family diversification.
- Barrio-Hernandez et al. 2023. https://doi.org/10.1038/s41586-023-06510-w — Scalable structural clustering; structural similarity needs evolutionary corroboration.
- Garg & Hochberg 2025. https://doi.org/10.1093/molbev/msaf124 — Empirical 3Di substitution model; proposed branch-rate application needs validation.
- Mutti et al. 2025. https://doi.org/10.1093/molbev/msaf149 — Benchmarks caution against replacing sequence phylogenomics by default.
- Szánthó et al. 2025. https://doi.org/10.1038/s41559-025-02851-z — Broad fungal sampling, close relatives, topology and calibration uncertainty. Its fungal boundary must be distinguished from database taxonomy.

These are methodological precedents, not an exhaustive systematic review. Access and update dates belong in source receipts for downloaded datasets.

- Galindo et al. 2021. Phylogenomics of a new fungal phylum reveals multiple waves of reductive evolution across Holomycota. https://doi.org/10.1038/s41467-021-25308-w — Published Sanchytriomycota genomes and proteomes in Figshare project 91439 permit filling an annotation gap.
- Grau-Bové et al. 2017. Dynamics of genomic innovation in the unicellular ancestry of animals. https://doi.org/10.7554/eLife.26036 — Genome-backed unicellular holozoan outgroup candidates and annotation bundles.

- Official AFDB v6 release notes: https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes — Current release differs from v4 bulk archives; unchanged coordinates may be relabeled. Model-version and sequence provenance must be retained.
- Official bulk archive description: https://github.com/google-deepmind/alphafold/blob/main/afdb/README.md — Taxid-sharded public bucket inventory. Observed suffix versions, not bucket name, determine selected archive version.
# Additional orthology method

- **OrthoFinder v3 (2026), improved phylogenetic orthology inference with enhanced accuracy and scalability**, Nature Methods, [DOI: 10.1038/s41592-026-03126-6](https://doi.org/10.1038/s41592-026-03126-6). Relevant to family inference and reconciliation at this project's scale. The [official advanced tutorial](https://orthofinder.github.io/OrthoFinder/tutorials/advanced-tutorial/) documents a diverse reference-core analysis followed by assignment of additional species and combined phylogenetic analysis. This is a computational decomposition of the complete cohort, not a separate pilot experiment. The official release selected for installation is [v3.1.5](https://github.com/OrthoFinder/OrthoFinder/releases/tag/v3.1.5), with its archive verified against the GitHub release API SHA256 digest. Orthology analyses have not yet been launched.

- Põlme et al. FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. Fungal Diversity 105, 1–16 (volume year 2020; online 2021). https://doi.org/10.1007/s13225-020-00466-2. Publisher genus-trait supplement supplies ranked candidate evidence. Genus matches are not species-validated traits, and the downloaded workbook is pinned separately from the historical article counts.
- Hess et al. 2018. Rapid Divergence of Genome Architectures Following the Origin of an Ectomycorrhizal Symbiosis in the Genus Amanita. https://doi.org/10.1093/molbev/msy179. Primary genomic comparison supplies focal species classifications and a within-genus ecological contrast; related symbiotic tips are not independent origins. Useful for taxonomy-sensitive trait curation and later gene-family analyses.
- Martin et al. 2008. The genome of Laccaria bicolor provides insights into mycorrhizal symbiosis. https://doi.org/10.1038/nature06556. Primary evidence for the focal species' ectomycorrhizal classification; it does not independently validate the phenotype of every subsequently sequenced isolate.

- van Kempen et al. 2023. [Fast and accurate protein structure search with Foldseek](https://doi.org/10.1038/s41587-023-01773-0). Native coordinate-derived 3Di encoding used in the executed benchmark. The specific build/source implementation was audited for invalid-state handling and six-residue feature context; structural search/encoding does not by itself validate an evolutionary substitution model.

- Garg & Hochberg model deposit, [Edmond DOI 10.17617/3.1MJJBH](https://doi.org/10.17617/3.1MJJBH), version 3.0. Executed use now includes the authors' coordinate-trained Q.3Di.AF model and ProstT5-trained Q.3Di.LLM sensitivity matrix; publisher checksums and exact files are recorded in metadata/3di_substitution_model_receipt.json. The prior proposed application has advanced to conditional point estimates, while adequacy and branch-uncertainty validation remain pending.

- Felsenstein 1985. [Confidence limits on phylogenies: an approach using the bootstrap](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x). Character resampling provides the reference for the paired site-resampling analysis. Its independence assumptions are not automatically satisfied by overlapping/spatially linked 3Di features; fixed-topology branch intervals here are not bootstrap clade support.

- NCBI. [GFF3 format](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/) and [The Genetic Codes](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/) (accessed 2026-09-13). Primary format references supporting documented table-1 defaults, explicit nonstandard codes, source-region code attributes and strand-aware partial-boundary interpretation. Applied in the marker annotation audit; annotation metadata does not independently prove biological completeness or sequence accuracy.

- Lofgren et al. 2021. [Comparative genomics reveals dynamic genome evolution in host specialist ectomycorrhizal fungi](https://doi.org/10.1111/nph.17160). Table 1 supplies species lifestyles and compiled Suillus host categories; uncertain and multiple hosts require explicit coding and independent transition reconstruction.
- Kohler et al. 2015. [Convergent losses of decay mechanisms and rapid turnover of symbiosis genes in mycorrhizal mutualists](https://doi.org/10.1038/ng.3223). Primary comparative genomics and root-expression evidence; Piloderma naming differs across later sources and remains flagged.
- Peter et al. 2016. [Ectomycorrhizal ecology is imprinted in the genome of the dominant symbiotic fungus Cenococcum geophilum](https://doi.org/10.1038/ncomms12662). Supplies species classifications and saprotrophic comparators; selected-isolate matching and current-panel transition independence remain unresolved.
- Martin et al. 2010. [Périgord black truffle genome uncovers evolutionary origins and mechanisms of symbiosis](https://doi.org/10.1038/nature08867). Primary Tuber melanosporum genome and ectomycorrhizal classification.

- Borneman et al. 2012 (online 2011). [The genome sequence of the wine yeast VIN7 reveals an allotriploid hybrid genome with Saccharomyces cerevisiae and Saccharomyces kudriavzevii origins](https://doi.org/10.1111/j.1567-1364.2011.00773.x). Primary selected-strain hybrid evidence; informs homeolog and species-tree sensitivity handling.
- Salazar et al. 2019. [Chromosome level assembly and comparative genome analysis confirm lager-brewing yeasts originated from a single hybridization](https://doi.org/10.1186/s12864-019-6263-3). Selected CBS 1483 genome study; hybrid ancestry requires explicit handling before conventional species-tree reconciliation.

- Tien et al. 2013. [Maximum Allowed Solvent Accessibilites of Residues in Proteins](https://doi.org/10.1371/journal.pone.0080635). Table 1 supplies theoretical ALLOWED-region reference areas and reproduces the older Miller scale. Terminal residues and ASA implementation differences require care; applied here as explicit normalization conventions with unclipped outputs.

- Sankoff 1975. [Minimal Mutation Trees of Sequences](https://doi.org/10.1137/0128004). Tree-based minimization of character-change cost. Here a unit-cost fixed-topology diagnostic is used, with missing tips unconstrained; it is not a fitted substitution-rate or ancestral-sequence uncertainty analysis.

- Naser-Khdour et al. 2019. [The Prevalence and Impact of Model Violations in Phylogenetic Analysis](https://doi.org/10.1093/gbe/evz193). Introduces maximum matched-pairs diagnostics used in IQ-TREE for substitution-model assumption screening. Applied without filtering; unavailable results and calibration limits remain explicit for short, masked 3Di alignments.

### Dating calibration provenance

Parham et al. (2012), *Best Practices for Justifying Fossil Calibrations*.
https://doi.org/10.1093/sysbio/syr107.
Specimen-level placement and geological age justification guide the calibration
review fields in dating-workflow.md. This methodological source does not supply
approved fungal calibrations for the project.

The existing Szánthó et al. (2025) archive was additionally audited for its four
mean chronograms and project name overlap; see dating-workflow.md and
metadata/published_chronogram_inventory_receipt.json. No published ages were
transferred to project nodes.

## Mixed-model computation

Bates D, Mächler M, Bolker B, Walker S (2015). Fitting Linear Mixed-Effects
Models Using lme4. *Journal of Statistical Software* 67(1):1–48.
[Primary article](https://www.jstatsoft.org/article/view/v067i01). Describes
profiled ML/REML computation and constrained covariance-parameter estimation.
Relevant to the proposed shared-control, family and phylogenetic mixed model;
it does not validate our particular covariance assumptions. Our specialized
evaluator is checked separately against direct dense calculations; optimization
and inferential calibration remain outstanding.

## Small-sample fixed-effect uncertainty

- Kenward & Roger (1997). [Small sample inference for fixed effects from restricted maximum likelihood](https://pubmed.ncbi.nlm.nih.gov/9333350/). Abstract reviewed; the original full text has not been reviewed here. Proposes adjusted fixed-effect covariance and an approximate F reference distribution. This is a candidate method, not evidence of coverage for the fungal matched model.
- Halekoh & Højsgaard (2014). [A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models – The R Package pbkrtest](https://www.jstatsoft.org/article/view/v059i09). Full author article downloaded and Appendix A inspected. Gives covariance adjustment and moment matching for Gaussian mixed models with independent, constant-variance residuals. Its known-kernel representation is relevant to our covariance model. Download provenance is in `metadata/mixed_model_uncertainty_literature_20260928.json`; project-specific feasibility and unresolved validation are in `matched-model-calibration-20260928.md`. This method has not been adopted for reported intervals.

## Ancestral insertion/deletion reconstruction

Ashkenazy H et al. (2012). FastML: a web server for probabilistic reconstruction
of ancestral sequences. *Nucleic Acids Research*40:W580–W584.
[Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3394241/)
([DOI](https://doi.org/10.1093/nar/gks498)). Describes separate binary coding
and reconstruction of multi-position gap characters alongside amino-acid
reconstruction, with ancestral uncertainty outputs. Relevant to avoiding an
implicit residue at every ancestral alignment column. It does not establish
that our alignments or inferred gap histories are correct. Local source
FastML3.11 was retrieved from the authors'
[source archive](https://fastml.evolseq.net/source/FastML.v3.11.tgz);
checksum and build plan are in `metadata/fastml_indel_build_plan_20260927.json`.
The public source page states a200-sequence limit, so local-tool capacity for
our622-protein family must be checked before scientific use. The wrapper's
default inclusion of terminal gaps requires explicit sensitivity analysis,
especially for extracted domains and incomplete annotations.

- Holmes 2017. [Historian: accurate reconstruction of ancestral sequences and evolutionary rates](https://doi.org/10.1093/bioinformatics/btw791). Candidate explicit indel/substitution framework. Pinned author implementation built and upstream tests passed; root, profile approximation and full input compatibility remain unqualified. See the [implementation assessment](historian-method-assessment-20260927.md).
- Westesson, Lunter, Paten and Holmes 2012. [Accurate reconstruction of insertion-deletion histories by statistical phylogenetics](https://doi.org/10.1371/journal.pone.0034572). Methodological predecessor cited by Historian. The implementation's nonoverlapping within-branch indel approximation must remain explicit; software validation is not replication of published benchmarks.


- Minh, Hahn and Lanfear 2020. [New Methods to Calculate Concordance Factors for Phylogenomic Datasets](https://doi.org/10.1093/molbev/msaa106). Gene concordance with variable taxon coverage complements clade support. Applied to all 250 marker trees and eight candidate species views; full independent readback and original-journal closure passed. Concordance remains conditional on gene-tree inference and does not identify causes of discordance. [Official IQ-TREE documentation](https://iqtree.github.io/doc/Concordance-Factor), accessed October 2 UTC, 2026; actual installed build and native output formats are pinned and validated separately.


- Zhang, Rabiee, Sayyari and Mirarab 2018. [ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees](https://doi.org/10.1186/s12859-018-2129-y). Full gene-tree species estimation with unresolved branches; informs the 30-run alignment/taxon/support sensitivity design. Published low-bootstrap contraction does not numerically validate our distinct SH-aLRT settings. Installed 5.7.8 jar/dependencies were checked against the official archive.
- Sayyari and Mirarab 2016. [Fast Coalescent-Based Computation of Local Branch Support from Quartet Frequencies](https://doi.org/10.1093/molbev/msw079). Method behind ASTRAL local support/coalescent internal lengths; locality and inferred-gene/model assumptions require qualification before treating support as calibrated topology uncertainty. Independent local/global quartet and numerical validation passed all 30 full candidates, followed by source/two-journal closure; biological qualification remains open.

- Sukumaran and Holder 2010. [DendroPy: a Python library for phylogenetic computing](https://doi.org/10.1093/bioinformatics/btq228). Tree parsing, manipulation and split-based comparisons underpin the independent raw-tree pruning/RF/bipartition readers. The installed implementation is checked using complete source trees and altered-output software contracts; software agreement does not qualify biological interpretation.
- Moreno, Holder and Sukumaran 2024. [DendroPy 5: a mature Python library for phylogenetic computing](https://doi.org/10.21105/joss.06943). Current major-version software reference. This project retains its installed 5.0.8 environment for the launched readers; no package upgrade is made during execution.

- ASTRAL author implementation, accessed October 2 UTC 2026: [WQInference.java](https://github.com/smirarab/ASTRAL/blob/master/main/phylonet/coalescent/WQInference.java) and [Posterior.java](https://github.com/smirarab/ASTRAL/blob/master/main/phylonet/coalescent/Posterior.java). Branch annotation retains available-gene effective N for resolved-evidence discrepancies at most0.001; otherwise it substitutes fractional resolved evidence. Explicit N enters fractions, posterior and MAP calculations. Upstream source inspection is supported by three actual installed5.7.8 jar threshold contracts, not treated alone as proof of build behavior. Full source/download provenance stays outside Git and is bound in the corrected audit plans.

## Full sequence-derived correspondence controls

- MAFFT author [manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html)
  and [acceptable-symbol documentation](https://mafft.cbrc.jp/alignment/software/anysymbol.html),
  accessed October 2 UTC 2026. Document `--auto` strategy selection, explicit
  amino-acid input and unknown-character scoring/preservation. The manual
  identifies its older coverage, so actual installed 7.525 commands and byte
  preservation are checked separately. These sources inform the full native
  sequence correspondence control, not a claim of correct project homology.
- FAMSA [author repository and usage](https://github.com/refresh-bio/FAMSA),
  accessed October 2 UTC 2026. Single-input MSA and explicit thread settings
  guide the second correspondence method. Installed 2.5.2-2598410 version/help,
  executable hashes and all six-order native software contracts are recorded
  separately. Cross-method agreement is a sensitivity diagnostic and does not
  validate homology or make sequence-derived predicted coordinates independent.
