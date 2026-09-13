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
