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
