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
