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

## Taxon identity and hybrid sensitivity

Name-based review flags 22 of the 526 working entries (21 incompletely identified fungal labels and one explicit Saccharomyces hybrid). Distinct taxon IDs must not be equated with established distinct species. Preserve current data and exploratory runs, but resolve names/assembly identity and assess exclusion of flagged taxa before final species-level inference. In particular, hybrid ancestry can violate a single bifurcating species-tree model. The hybrid is involved in several large exploratory global RMSDs; this is a review priority, not evidence that hybridization caused structural change. See `metadata/taxon_label_review.tsv` and reproduce with `python scripts/audit_taxon_labels.py`.

## Paired structural branch inference

Use identical observed cells and sequence-based residue correspondence for AA and native coordinate-derived 3Di fits. Preserve invariant sites and missing-data reasons. Fit structural branch lengths on the sequence marker topology using the published AF empirical model, with frequency/training-source sensitivities. Preserve unrooted splits and full taxon-set identity; do not compare edges from different taxon subsets as though they are the same branch. Treat near-zero lengths and rare states as uncertainty/model concerns, and withhold branch ratios and acceleration rankings until those concerns are addressed.

## Dependence and conditional resampling

Resample identical columns in AA and 3Di alignments across all taxa. Use 200 draws each at block lengths 1, 10 and 30 as a conditional sensitivity analysis, retaining invariant sites and recording unestimable draws rather than silently redrawing. Do not interpret percentile intervals as overall calibrated confidence: the audited inputs contain many spatially linked feature pairs farther apart than the chosen blocks. Preserve paired sampling covariance for later error-aware coupling analysis; it is not an evolutionary correlation.

## Source-stratified structural snapshots

Use an explicitly selected provider and prediction tool for the primary evolutionary mapping; rank candidates only within that source and leave unavailable sources missing. Retain all alternative models for same-sequence source comparisons. Higher pLDDT across pipelines is not evidence of greater accuracy and must not silently switch the production source. The first available GDM/ColabFold pair is a control with limited scope, not a calibration of pipeline effects across fungi.

## Cloud availability and prediction strategy

The user confirmed that no existing Google Cloud account is available and asked whether ESMFold would be faster. Cloud authentication is therefore not an available acceleration route for the previously attempted proteome archives. Continue the already authorized individual-model downloads and local ESMFold execution; no account creation or paid cloud provisioning is planned.

At this check, the active ESMFold run had approximately 3,400 completed per-sequence prediction receipts, with median measured inference approximately 1.69 seconds and mean 2.43 seconds for the completed subset. The current production queue is limited to canonical proteins of at most 512 residues; this is a local execution limit, not an ESMFold hard limit. Completed short-sequence timings do not forecast long proteins or the entire multi-million-protein atlas. Individual AlphaFold acquisition and PAE acquisition are also progressing, so replacing already available models with new predictions would duplicate work. Continue complementary retrieval and prediction, with explicit predictor strata and same-sequence controls to avoid confounding clade effects with prediction methods. EBI also documents public FTP bulk subsets; the earlier authentication failure concerned the specific Google-hosted archives attempted, not all AlphaFold access.
