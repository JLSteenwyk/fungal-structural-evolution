# Metadata

`assembly_candidates.tsv`: all catalog records retained with source label; not a selected sample. Assembly accession includes version; species_taxid identifies NCBI species grouping; taxid can identify a strain. Paired GCA/GCF assemblies are not independent species.

`sampling_manifest.tsv`: selected/review taxa only. Fields: taxon_id (stable project ID), species_name, species_taxid, study_role (ingroup/outgroup/boundary_review), lineage, assembly_accession, annotation_source_version, proteome_url, cds_url, genome_url, assembly_level, busco_complete_pct, contamination_status, ecology, ecology_source, structure_coverage, inclusion_reason, status. Empty values mean unknown, never zero or absence. No rows yet because selection is not complete.

`source_receipts.json`: source URL, retrieval time, local ignored path, SHA256 and byte size. Catalog hashes pin snapshots; remote content may change. Preserve local snapshots or archive externally before publication.
