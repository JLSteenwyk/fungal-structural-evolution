# Project data dictionary

This index defines the primary project-level tables and immutable receipts. It
does not replace stage-specific dictionaries, which remain the authoritative
schema for their corresponding analysis outputs. Large data, coordinates, PAE
matrices, and structure-search databases are held outside Git; their locations,
checksums, and retrieval provenance are recorded by the receipts linked below.

## Cohort manifests

`metadata/analysis_manifest.tsv` is the frozen 526-taxon analysis cohort: 501
fungal ingroup entries and 25 non-fungal outgroups. `metadata/sampling_manifest.tsv`
is the 527-entry candidate cohort and additionally retains *Saccharomyces
jurei*, which is excluded from analysis because no usable annotated proteome
was acquired. Do not substitute one manifest for the other without recording
the cohort change.

Both use these columns:

| Field | Meaning |
| --- | --- |
| `taxon_id` | Stable project taxon identifier; use this key for joins. |
| `species_name` | Literal name retained in the project record. |
| `species_taxid` | Source taxonomy identifier, where available. |
| `study_role` | `ingroup` or `outgroup`; this is a study designation, not a taxonomic claim. |
| `lineage` | Semicolon-delimited retained lineage label. |
| `assembly_accession` | Assembly accession or externally sourced assembly identifier. |
| `annotation_source_version` | Source and version/date for the annotation used. |
| `proteome_url`, `cds_url`, `genome_url` | Recorded retrieval locations for proteins, coding sequence, and assembly. |
| `assembly_level` | Assembly-level descriptor supplied by the source. |
| `busco_complete_pct` | Recorded BUSCO completeness percentage; blank means unavailable. |
| `contamination_status` | Recorded contamination-review disposition. |
| `ecology`, `ecology_source` | Trait assignment and its source; `unknown` or blank remains explicit missingness. |
| `structure_coverage` | Initial structure-coverage disposition; it is not a downstream eligibility call. |
| `inclusion_reason` | Reason the taxon was retained in the candidate design. |
| `status` | Current cohort/QC disposition. |

The membership, 25-outgroup balance, source hashes, and scope limitations are
independently checked in
`metadata/outgroup_selection_audit_20261006_v1.json`. See
[outgroup selection](outgroup-selection-20261006.md) for the selection
rationale.

## Orthology input inventory

`metadata/orthology_input_manifest.tsv` records the protein FASTA supplied to
the scalable orthology workflow. Its fields are `taxon_id`, `species_name`,
`study_role`, `lineage`, `stage`, `proteins`, `source_sha256`, `input_path`,
and `reason`. `source_sha256` binds the exact input file; `stage` describes
the workflow partition and does not denote an evolutionary classification.
See [the orthology workflow](orthology-workflow.md) for family and
reconciliation semantics.

## Atlas-availability and uncertainty records

`metadata/full_prediction_atlas_union_taxon_coverage_20261005_v1.tsv` gives
taxon-level availability counts. `representative_proteins` is the denominator;
`afdb_only`, `esmfold_only`, `both`, and `neither` are mutually exclusive
source-availability states. Availability is not a confidence filter, a
homology call, or structural-analysis eligibility.

The PAE recovery preparation receipt,
`metadata/full_atlas_missing_pae_recovery_preparation_20261021_v1.json`, binds
the interrupted source queue and recovery queue with hashes and records
`expected_models`, retained disposition counts, source-state census, and
script identity. Its completion receipt is intentionally absent while recovery
is live. The queued readback, explicit-failure retry, retry readback, and
lineage composition each produce separate immutable receipts; the full atlas
cannot enter clustering until the [gate report](../scripts/report_full_atlas_gate_status_v2.py)
is ready. See [the PAE workflow](full-atlas-pae-20261005.md).

## Phylogeny plan and model outputs

`metadata/species_taxon_pmsf_plan_20261001.json` is the pinned plan for the
16 crossed native PMSF refits. `inputs`, `output`, `executable`, `seed`,
`cached_guides`, `jobs`, `pins`, `resources`, and `scope` define its frozen
inputs, model grid, computational limits, and interpretation limits. A run
receipt records computational completion only after all required outputs and
audits pass; it is not itself selection of a preferred species tree. See the
[phylogenetic workflow](phylogenetic-workflow.md).

## General rules

- Join tables by their documented stable identifiers and retain unmatched rows
  as explicit dispositions.
- Treat blank fields, `unknown`, and review states as missing or unresolved;
  do not silently recode them as absence or a biological result.
- Use source hashes and immutable receipts to reproduce an executed stage.
  Do not use a changing runtime state file as a final analysis record.
- Structural clusters, model availability, and prediction confidence do not
  define orthogroups, experimentally established function, or selection.
