# Metadata

`assembly_candidates.tsv`: all catalog records retained with source label; not a selected sample. Assembly accession includes version; species_taxid identifies NCBI species grouping; taxid can identify a strain. Paired GCA/GCF assemblies are not independent species.

`sampling_manifest.tsv`: selected/review taxa only. Fields: taxon_id (stable project ID), species_name, species_taxid, study_role (ingroup/outgroup/boundary_review), lineage, assembly_accession, annotation_source_version, proteome_url, cds_url, genome_url, assembly_level, busco_complete_pct, contamination_status, ecology, ecology_source, structure_coverage, inclusion_reason, status. Empty values mean unknown, never zero or absence. No rows yet because selection is not complete.

`source_receipts.json`: source URL, retrieval time, local ignored path, SHA256 and byte size. Catalog hashes pin snapshots; remote content may change. Preserve local snapshots or archive externally before publication.

`assembly_candidates_taxonomy.tsv`: NCBI ranked lineage joined by species_taxid. Missing IDs remain missing; taxonomy updates and merged IDs require explicit follow-up.

`taxonomy_coverage.json`: species coverage per NCBI phylum, including unclassified and missing-lineage cases. NCBI Fungi membership does not resolve competing biological circumscriptions.

`species_assembly_candidates.tsv`: one provisionally preferred annotated latest assembly per species. URL templates are not yet availability-verified. This is not the final sample.

`timetree_published_taxa.tsv`: 153 taxa and original data sources from Szánthó et al. Supplementary Tables 1–2. Published roles and taxonomy do not assign project roles; some published sources are transcriptomes rather than genomes. `timetree_source.json` preserves publisher file metadata and license. The downloaded archive's publisher MD5 was verified.

Use `make restore-sources` to verify cached snapshots or restore them if the remote content still matches recorded SHA256. If mutable NCBI content has changed, the command fails rather than claiming exact reproduction. `make taxonomy`, `make published-taxa`, and `make species-candidates` rebuild derived tables. openpyxl is needed for the published workbook import.

`resolved_species_taxid` and `taxonomy_id_status` preserve current/merged/deleted/unresolved distinctions. Original assembly species IDs are retained. Deleted IDs are not automatically substituted.

`outgroup_assembly_candidates.tsv`: targeted name-based discovery from NCBI protozoan catalogs; requires taxonomic and annotation review and is not the final outgroup sample. Rebuild with `python scripts/inventory_outgroup_sources.py` after restoring sources.

`proteome_availability.tsv`: generated when the resumable `scripts/check_proteome_availability.py` run completes. HEAD response and check time only; no claim that FASTA content or completeness was validated. Raw checkpoint records are ignored data, keyed by URL; malformed historical URLs are not used for corrected candidate URLs.

`fungal_sampling_draft.tsv` and its summary describe 500 unique fungal species selected for full-scale QC, not the final frozen sample. The main sampling_manifest remains unpopulated until review. Selection is reproducible from config/sampling.json, with unknown ecology retained explicitly.

`proteome_download_receipts.json`: completed download/QC records (can be partial during execution), publisher MD5, SHA256, source URL, local path and protein statistics. `proteome_qc_snapshot.json`: aggregate of those records, not a liveness indicator. Download via `python scripts/download_proteomes.py`; derive progress via `python scripts/summarize_downloads.py`. FASTA validity is not biological completeness.

`fungal_sampling_expanded_draft.tsv`: initial 500 plus two Sanchytriomycota species using published annotations. Total 502, consistent with approximate target. Some original candidates require replacements after availability/QC; this remains provisional. Rebuild via scripts/add_sanchytriomycota.py after external retrieval.

`figshare_*.json`: pinned publisher versions, file lists, licenses and checksums for external genome bundles. `external_genome_receipts.json`: successfully verified files so far; partial while acquisition runs. `scripts/fetch_external_genomes.py` resumes from recorded publisher metadata and verifies cached content. Genome/proteome presence does not prove assembly identity across repositories.

`qc_input_receipts.json`: raw proteomes normalized for preliminary BUSCO completeness. Terminal `*` and `.` markers removed; sequences with internal stops or other noncanonical characters excluded and counted. Source and output hashes preserved; no isoform collapse yet. This is separate from final orthology inputs.

`busco_dataset_receipt.json`: dataset configuration and SHA256 for every file in the initial eukaryota_odb12.2 marker set. `environments/busco-linux-64.lock.txt` pins installed packages. Invoke BUSCO via conda run with .cache/envs/busco; directly invoking its script under the base Python does not activate the environment.

`outgroup_sampling_draft.tsv`: exactly 25 genome-backed taxa; selection logic and citations in config/outgroups.json. External `assembly_accession` values prefixed `figshare:` identify the actual publisher assembly artifact, not an asserted NCBI match. Missing species taxids are unresolved identifiers, not absent organisms. `outgroup_ncbi_download_manifest_download_receipts.json` records all 20 validated NCBI outgroup downloads.

`sampling_manifest.tsv` is now the combined provisional 527-taxon manifest, superseding earlier notes that it was empty. It includes the unresolved S. jurei record; only 526 taxa currently have acquired protein data. Rebuild with scripts/build_sampling_manifest.py. Final annotation QC, ecological metadata and structure coverage remain pending.

`busco_eukaryota_qc.tsv`: successful raw-proteome BUSCO runs, counts checked for internal consistency, source summary hashes retained. Rebuild with scripts/summarize_busco.py. It does not represent lineage-specific completeness or corrected isoform counts.

`analysis_manifest.tsv`: 501 available fungal species plus 25 outgroups; candidate exclusion documented separately. Still subject to biological QC and trait review.

`alphafold_bulk_coverage.tsv`: archived bulk availability at queried species and selected assembly taxids. `alphafold_bulk_objects.json` retains object generations, update dates, sizes and publisher hashes. `alphafold_selected_archives.json` retains highest-version-per-shard objects. No exact protein-sequence coverage is implied; no structure archive downloaded in this inventory. Reproduce with inventory_alphafold_archives.py then select_alphafold_archives.py. Query cache resides in ignored data/raw.

`data/raw/uniprot_matches.jsonl` (ignored): resumable per-taxon exact matching records. Per-taxon reference tables, release headers, checksums and protein-to-UniProt mappings reside under data/uniprot. Script match_uniprot_structures.py requires the combined analysis and QC input manifests. AlphaFold crossrefs are candidate availability evidence only; sequence equality to the actual downloaded model remains required.

Legacy archive retrieval code is preserved but was stopped after verified access denial. It supports generation-pinned range resume and publisher MD5 verification if authorized object access becomes available. Do not mistake listed object metadata for a downloaded structure.

`structure_reuse_snapshot.json` and `verified_structure_receipts.json`: partial current-AFDB retrieval evidence, derived by summarize_structure_reuse.py. Full raw API files and CIFs live in data/structures/afdb; protein-level mappings are in input_links.tsv. Multiple proteins may share one model only on exact sequence identity. Confidence summaries are measurements, not blanket acceptance of structural accuracy. No PAE-based domain orientation analysis has occurred.

`taxon_label_review.tsv` records explicit hybrid, incompletely identified and duplicate name flags from the working manifest. `flags` is a semicolon-separated list; `review_status` remains pending until taxonomy is independently assessed. This screen does not establish that unflagged entries are distinct accepted species. Its receipt pins the source manifest and output checksum.

`mafft_alignment_statistics.tsv` and its audit receipt summarize all 125 completed full-protein alignments; `mafft_species_matrix_receipt.json` pins the 526-taxon alternative concatenation. `alignment_correspondence.tsv` reports shared retained residue counts, each method's cross-taxon residue-pair edges, shared edges, edge Jaccard (intersection/union), and directional recovery fractions. These are conditional on residues retained by both 50%-occupancy masks, not accuracy estimates. Its receipt pins individual alignment receipts, totals and detailed result checksums. Per-column diagnostics remain outside Git under the recorded result snapshot.

`pfam_release_receipt.json` pins Pfam 38.2 resource URLs, publisher MD5s, compressed/decompressed SHA256 hashes, version text and all 64 HMM chunks. `marker_domain_input_receipt.json` records exact-sequence deduplication and hashes of the full marker protein links. `marker_domain_search_config.json` records the active HMMER binary hash/version, worker settings, threshold rule and source receipts. These preparation/configuration receipts do not establish completed domain annotations; completion requires all search chunks and the separately validated annotation summary.

`assembly_quality_metrics.tsv` records verified NCBI statistics, with `all_` and `primary_` prefixes identifying different assembly scopes; missing values remain blank. External taxa are explicitly directed to the separate table. `external_assembly_quality.tsv` reports direct FASTA-record statistics for all seven external genomes; its record N50 must not be silently merged with NCBI contig N50. The respective receipts pin sources, outputs and scripts. `assembly_quality_figure_receipt.json` links both to the full 526-taxon BUSCO comparison; raw reports and individual provenance remain at its recorded data paths.

`ecology_candidates.tsv` contains genus-ranked evidence only; `species_verified` and `eligible_for_confirmatory_tests` are false throughout. Its receipt pins the publisher workbook, exact-match policy, duplicate names and source-row counts. `species_ecology_evidence.tsv` contains separately curated primary-literature statements from `config/ecology_species_evidence.json`, with species/variety scope, source locators, selected assembly and conflicts with genus projections. Published species classifications are distinct from exact-isolate experimental verification and from eligibility for a replicated-transition test. Neither table automatically overwrites the main manifest's ecology fields.

`full_domain_input_receipt.json` records complete representative-proteome inputs for domain search, including every taxon's source checksum/count, exact unique-sequence totals, marker-query reuse and additional FASTA/link hashes. `full_domain_search_estimate.json` provides a provisional runtime projection with the observed chunk receipts and input/configuration hashes. Neither receipt proves full-proteome domain annotation is complete. Large search inputs and all protein/taxon links remain under `data/domains/full-inputs-v1/`.

- `marker_domain_search_receipt.json` and `marker_domain_annotation_receipt.json`: completed all-profile marker searches and validated raw hits; overlaps retained, not resolved architectures.
- `esmfold_checkpoint_receipt.json`: pinned publisher revision, weight checksum and model-file identities; local verified cache reuse.
- `marker_prediction_input_summary.json`: frozen queue counts and artifact hashes; full per-taxon match hashes remain in the ignored input receipt. Inventory/retrieval pending is distinguished from absence of a reuse candidate.

- `local_prediction_chunk1.tsv` and `local_prediction_chunk1_receipt.json`: independently read-back first 128 predictions, measured resources and confidence, with coverage and remaining-work counts.
- `esmfold_execution_audit.json`: checkpoint-loading exceptions, rejecting contact-head hook, repeated inference and actual ESM Python source hashes.
- `orthology_core_receipt.json`: successful 64-taxon computational core; not evidence of full-cohort orthology completion.

- `orthofinder_numpy_install_report.json` and `orthofinder_numpy_compatibility.json`: NumPy 2.2.6 wheel provenance and successful round-trip through the assignment function that failed under 2.5.3; full install report outside Git.

- `structural_domain_comparison_receipt.json`, `structural_domain_figure_receipt.json` and `structural_domain_review_order.tsv`: 738 conserved-domain comparisons, source hashes, plot provenance and complete manual review ordering. Refitting improvement is not evidence of biological motion.
- `domain_placement_receipt.json` and `domain_placement_confidence.tsv`: 375 domain-pair/taxon-pair combinations at three directional PAE thresholds; absent confident residue-pair measurements remain missing, not zero.

- `foldseek_native_config.json`: pinned binary/source identities and native coordinate-extraction commands, including recorded header-link recovery.
- `native_3di_reproduction_receipt.json`: exact semantic agreement between the original native extraction and the complete scripted reproduction for all 422 models.
- `3di_feature_audit_receipt.json` and `3di_model_summary.tsv`: native feature/partner validity and focal/six-residue confidence masks, with hashes of the per-model NPZ artifacts.
- `3di_geometry_benchmark_receipt.json`, `3di_benchmark_figure_receipt.json` and `3di_confidence_summary.tsv`: matched-site geometry/state benchmark provenance and counts under three confidence regimes; not branch estimates.
