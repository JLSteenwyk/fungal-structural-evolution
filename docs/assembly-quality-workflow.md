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
