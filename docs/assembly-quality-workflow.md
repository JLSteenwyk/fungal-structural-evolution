# Assembly contiguity and metadata quality

Assembly statistics supplement marker recovery and annotation reconciliation. They do not establish absence of contamination, correct gene models or independent species identity.

## NCBI assembly reports

`retrieve_assembly_statistics.py` retrieves the `_assembly_stats.txt` report belonging to each of the 519 selected NCBI assembly versions. The exact accession must appear in the report's GenBank/RefSeq header. Publisher MD5 values and local SHA256 hashes are verified, and raw reports plus per-assembly receipts remain under ignored `data/assembly_statistics/`. Completed cached files are hash-checked and reparsed. Two download workers are used.

The parser keeps statistics for `all/all/all/all` separate from `Primary Assembly/all/all/all`. No component, organelle or primary-assembly value is silently substituted for a whole-assembly statistic. It rejects duplicate scope/statistic keys, invalid numbers, component lengths exceeding totals and invalid GC percentages. Absent fields remain missing; they do not become zeros. Assembly type, representation, sequencing method, reported coverage, BioSample and BioProject are retained as deposited metadata. The reported assembly type does not independently establish biological ploidy.

All 519 reports passed validation and total 3,834,008 bytes. Whole-assembly contig N50 is present for 474 reports (454 fungal entries and all 20 NCBI outgroups); 45 fungal reports omit that field. Primary-assembly scaffold N50 is present for 392 reports. Metric definitions follow the [NCBI assembly report documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/).

## Seven external genomes

`measure_external_assemblies.py` matches each external taxon's genome URL to a unique, checksum-verified deposited FASTA in the existing source receipts. It checks unique nonempty sequence IDs and IUPAC DNA characters, then computes record counts, total length, record N50/L50, longest record, N content, other ambiguity and GC percentage among A/C/G/T bases. All seven genomes passed.

These are **FASTA-record** statistics. Records are not split at gaps or assumed to be equivalent to NCBI contigs; organelles and all other deposited sequences remain included. Pirum gemmata and Abeoforma whisleri have record N50s of 1,637 and 2,170 bp, respectively, alongside 5.6% complete broad BUSCO-marker recovery each. Their fragmentation and sparse marker recovery motivate taxon-exclusion and gene/domain-detection sensitivities; they do not prove that any particular inferred absence is artifactual or evolutionary.

## Reproduction and interpretation

```bash
python scripts/retrieve_assembly_statistics.py
python scripts/measure_external_assemblies.py
python scripts/plot_assembly_quality.py --output results/qc/assembly-quality-v1
```

Use a new output directory for each figure snapshot. The figure verifies the two assembly-metric tables against their receipts and checks every BUSCO table value against its hashed source summary. Separate panels show NCBI whole-assembly contig N50 and external FASTA-record N50; missing N50 values are omitted from the scatter, not imputed. The joined table retains all 526 taxa and separate N50 columns. No regression, causal estimate or significance test is applied.

Summary tables and receipts are versioned under `metadata/`; raw reports and detailed figure outputs remain outside Git. Plot dependencies use `environments/structural-comparisons.yml`; acquisition uses Python's standard library and external FASTA measurement additionally uses Biopython. Tests cover assembly-version rejection, scope separation, missing-versus-zero values, duplicate metrics, invalid lengths and the exact N50 boundary. Contamination assessment, lineage-specific completeness, annotation-error review and taxonomic identity resolution remain outstanding.

## Selected-assembly FCS report inventory completed

Retrieved NCBI Foreign Contamination Screen reports from the exact assembly-version
directories for all 519 selected NCBI genomes. All 526 manifest entries remain in
the summary, including seven external genomes with no NCBI report requested.
The [NCBI FCS documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/contamination/fcs-contamination/)
describes these reports and distinguishes a header-only report with no flagged
sequences from missing evidence. The [publisher report specification](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_fcs_summary_and_details.txt)
defines one-based inclusive intervals and separate action labels; it is archived
with a checksum under `results/qc/fcs-source-documentation-v1`.

Validation found:

- 472 reports matching their publisher MD5 checksum.
- 46 reports absent from the publisher checksum listing, with identical contents
  verified by a second independent HTTPS request. Their weaker checksum
  provenance remains separately labeled.
- One unresolved publisher MD5 disagreement for F104408, GCA_013423385.1. Repeating
  both downloads reproduced the disagreement and identical report contents; that
  assembly's contamination evidence remains unresolved.
- Seven external genomes needing a separate contamination assessment.

Among the 518 usable NCBI reports, 86 contain EXCLUDE, FIX or TRIM regions. The
complete action-preserving table contains 4,039 region records. Counts are
publisher observations, not newly inferred biological contamination or confirmed
misannotations. No sequences or taxa were removed. A report without flagged rows
is not proof of absence of contamination, and FCS database/run versions remain
explicit in each report header. These historical FCS runs are not rescreens with
a current database.

The audit reads all original report rows back into the original region table and
independently recomputes interval unions with an event sweep. Overlapping spans
are counted once per sequence for each action group. EXCLUDE/FIX/TRIM are grouped
for review; REVIEW_RARE, REVIEW, INFO and organelle labels remain distinguishable.
Three parser tests passed, and 500 deterministic interval examples matched both
union implementations against explicit per-base sets. The raw report headers,
region fields, retrieval times, URL versions, publisher listings and checksums
are retained. Exact assembly URLs plus checksum evidence bind the files to the
requested versions; the FCS header itself does not supply an assembly accession.

```bash
python scripts/retrieve_selected_fcs_reports.py --manifest metadata/analysis_manifest.tsv --output results/qc/selected-assembly-fcs-v1
python scripts/audit_selected_fcs_inventory.py --inventory results/qc/selected-assembly-fcs-v1 --manifest metadata/analysis_manifest.tsv --output results/qc/selected-assembly-fcs-audit-v1
```

Use new output paths for reruns. The acquisition used two download workers;
the audit used one worker and repeated only missing-checksum reports. Large data
and per-assembly receipts remain outside Git. Compact receipts and the complete
526-entry summary are versioned as `metadata/selected_assembly_fcs_*`.
Mapping flagged genomic intervals onto CDSs, retained genes and structural
markers is the next step. It will establish which downstream comparisons need
contamination-related sensitivity analysis, without treating every assembly-level
flag as affecting every protein from that assembly.
