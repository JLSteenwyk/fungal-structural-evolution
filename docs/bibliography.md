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
