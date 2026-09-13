# Decisions and unresolved questions

- User-directed scope: approximately 500 unique fungal species and 25 additional non-fungal outgroups, no separate pilot.
- Repository name: fungal-structural-evolution. User confirmed public repository JLSteenwyk/fungal-structural-evolution; origin configured and remote verified empty before initial push.
- NCBI RefSeq and GenBank fungal assembly catalogs are initial discovery sources, not sufficient coverage by themselves. Augment sparse lineages with published datasets and repositories.
- Fungal circumscription is unresolved. Szánthó et al. (2025) use a narrower osmotrophic Fungi definition, treating Aphelida, Rozellida and Microsporidia as close relatives. Record NCBI classification separately from study role; do not silently count boundary taxa as non-fungal outgroups. Assess alternative definitions before freezing roles.
- Shared host: do not infer exclusive resource allocation from physical capacity. Use conservative metadata jobs initially; no GPU jobs launched.
- Existing `gh` executable is not the GitHub CLI (it rejects `auth status`); a GitHub configuration file exists. Do not expose credentials. Establish an appropriate client before remote creation.

Open: outgroup availability and boundary definition; trait sources/replicate transitions; predictor/database availability; actual reusable structure coverage; sustained shared-machine allocation.

- Source audit: Szánthó et al. Supplementary Table 1 assigns the Rozella allomycis JGI URL to Parvularia atlantis as well. Retain the original source row but do not download that URL as Parvularia. NCBI independently lists Parvularia assembly GCA_943704415.1, annotation availability pending.
- All 35 originally missing taxonomy IDs are in NCBI delnodes.dmp; none are merged IDs. Keep their records and flag deleted; do not infer replacements from similar names.
- Availability validation exposed trailing slashes in NCBI ftp_path. The initial malformed URL checks are invalid for taxon exclusion. Corrected URL construction strips trailing slashes; check cache is keyed by full URL, so corrected URLs are checked independently. HEAD availability remains distinct from validated FASTA content.

- Initial 500-genome QC draft uses the broad NCBI Fungi circumscription, including Microsporidia, Rozellomycota and Aphelidiomycota as ingroup. These cannot also count toward the 25 non-fungal outgroups. A narrower osmotrophic clade can be tested with explicit re-rooting/pruning and separate counts.
- Draft selection prioritizes unrepresented orders/families/genera and up to 12 annotated species from smaller phyla, plus focal genus sampling. This is a transparent taxonomic proxy, not branch-length-based phylogenetic optimization. Focal genera provide potential contrasts, not verified ecological transitions. Sanchytriomycota requires an external annotation source and is a known coverage gap.
- Full-scale proteome QC began directly for the 500-species draft. QC failures may require replacement; final sampling is not frozen.
