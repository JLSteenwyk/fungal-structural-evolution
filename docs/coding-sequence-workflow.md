# Coding-sequence acquisition and validation

The selection analyses require coding sequences tied to the same assembly, annotation and protein identifiers as the structural and orthology inputs. Acquisition has started for all **519 NCBI-backed taxa** in the 526-taxon working manifest. Every CDS URL must have the same assembly directory and prefix as its recorded protein FASTA. The prelaunch plan permits two download/validation workers, 50 GB additional disk and approximately 0.5–8 hours on the existing host; these are planning allowances, not measured completion costs. No paid resources are used.

```bash
python scripts/download_coding_sequences.py
python scripts/extract_creolimax_cds.py --output results/cds/creolimax-extraction-v1
python scripts/inventory_external_cds.py --creolimax results/cds/creolimax-extraction-v1
python -m unittest discover -s tests -p test_cds_fasta_validation.py
```

For each NCBI assembly, the downloader archives the publisher checksum listing, requires exactly one checksum for the CDS file, checks MD5 before accepting the download and records local SHA256. It checks gzip integrity, unique FASTA record identifiers, nonempty sequences and the IUPAC DNA alphabet. Counts of sequences without a protein identifier, non-triplet lengths and ambiguous bases remain explicit. Those properties are flags for later interpretation, not automatic deletion criteria. Three tests confirm that partial/unlinked records remain flagged, repeated IDs are rejected and non-DNA input fails validation.

The downloader holds a single-writer lock and resumes validated sources only after checking their sequence and checksum-listing hashes against the cached receipt and exact taxon/assembly/URL identity. SIGTERM or SIGINT stops new submissions and drains the two active downloads. Download errors and incomplete taxa remain in the event log and final status. Changed cached sources or publisher checksum listings require review rather than overwriting previously validated evidence. A completed download does not establish translation correctness.

Files remain under `data/cds`; each taxon's result is appended to `data/raw/cds_downloads.jsonl`. The completed batch summary is `metadata/cds_download_receipts.json`: all 519 NCBI taxa validated, with no pending/error taxa and the seven external-source taxa listed separately. Files occupy 2,768,796,083 compressed bytes. Independent receipt readback checks all CDS SHA256 and publisher MD5 values, archived checksum-listing hashes and exact manifest assembly/URL coverage; evidence is in `metadata/cds_acquisition_readback.json`. Execution log: `logs/cds_downloads_v1.log`. The launch plan is `metadata/cds_acquisition_resource_plan.json`.

## External-source taxa

The external-source inventory checks actual CDS file hashes, unique sequence IDs and the DNA alphabet against the existing published/extraction receipts. It now identifies available CDS sources for all seven taxa; translation-verified subsets and boundary-aware complete-codon projections remain distinct:

| Taxon | Available CDS records | Current status |
| --- | ---: | --- |
| Amoeboradix gromovi | 7,220 | Previously extracted with exact protein translation |
| Sanchytrium tribonematis | 9,367 | Previously extracted with exact protein translation; one source protein remains out of assembly bounds |
| Corallochytrium limacisporum | 7,535 | Genome-linked complete-codon spans reproduce every protein; partial boundaries recorded |
| Chromosphaera perkinsii | 12,463 | Genome-linked complete-codon spans reproduce every protein; partial boundaries recorded |
| Pirum gemmata | 21,835 | Genome-linked complete-codon spans reproduce every protein; partial boundaries recorded |
| Abeoforma whisleri | 17,283 | Genome-linked complete-codon spans reproduce every protein; partial boundaries recorded |
| Creolimax fragrantissima | 8,558 | Extracted with exact protein translation; 136 other proteins remain explicit exceptions |

Sources and artifact hashes are in `metadata/external_cds_source_inventory.json`. The two sanchytrid CDS sets inherit the detailed coordinate, strand and translation audit described in the gene/isoform workflow; their provisional ORF status and the Sanchytrium exception remain explicit. Raw external source files are preserved. The available published CDS counts do not prove one-to-one protein mappings or complete gene annotations.

## Requirements before selection tests

Next, link CDS records to exact versioned protein identifiers and gene-representative decisions. Validate translations under the documented genetic code, with terminal-stop handling, partial CDSs, ambiguous residues, frameshifts and annotation-specific exceptions explicit. Creolimax extraction has passed the strand, phase and exact-translation checks for its verified subset; its exceptions require further source review before inclusion. Never force translation agreement by trimming unexplained residues or picking a different code solely to improve the match.

Codon alignments must preserve validated amino-acid correspondence and gene-copy identity. Test suitability at the family/clade scale, including alignment quality and synonymous divergence, before choosing selection models. Full-genome downloads or raw sequence availability do not establish suitability for deep fungal comparisons. Selection tests, multiple-testing correction and mapping supported sites to structures remain pending; structural acceleration alone is not evidence of positive selection.

## Completed Creolimax CDS extraction audit

`extract_creolimax_cds.py` checks deposited genome, GTF and protein-source hashes, as well as the normalized protein input. It uses exact transcript IDs for 8,644 proteins and explicit gene-ID links for 50 assembler products. Gene-ID fallback is accepted only when one transcript is linked. All 8,694 input proteins have an audit row; ambiguous mappings remain exceptions.

CDS blocks are ordered in transcription direction, reverse-complemented on the minus strand, checked for overlapping intervals and assembly bounds, and joined with their split-codon bases intact. Internal GTF phase describes codon continuity; those bases must not be trimmed from every exon. This audit requires initial phase zero, valid internal phases and a total length divisible by three. Translation under table 1 must match the complete normalized protein after removing at most one terminal stop from the translation. A terminal stop remains in the nucleotide output and is flagged for later codon processing.

The completed results contain **8,558 exact translations**. The other **136 proteins** comprise 71 translation mismatches, 60 nonzero/missing initial phases, four non-triplet CDS lengths and one ambiguous gene-to-transcript mapping. They are not repaired by arbitrary trimming or alternate-code selection. Two annotation transcripts lack a selected protein association under these rules and remain listed separately.

```bash
python scripts/extract_creolimax_cds.py --output results/cds/creolimax-extraction-v1
python scripts/inventory_external_cds.py --creolimax results/cds/creolimax-extraction-v1
python -m unittest discover -s tests -p test_creolimax_cds.py
```

Four focused tests cover a codon split across exons, negative-strand ordering, inconsistent phase rejection and out-of-bounds rejection. Independent output readback reproduced all 8,558 translations and verified artifact hashes and the full 8,694-protein status partition. The verified FASTA, ordered annotation segments and unused-transcript list remain in `results/cds/creolimax-extraction-v1`; the receipt and complete exception rows are versioned as `metadata/creolimax_cds_extraction_receipt.json` and `metadata/creolimax_cds_exceptions.tsv`. This makes a verified CDS subset available for the seventh external taxon; codon-alignment and selection eligibility still require family-specific assessment.

## Published outgroup CDSs: strict and boundary-aware validation

All 59,116 published CDS identifiers map exactly to normalized proteins across Corallochytrium, Chromosphaera, Pirum and Abeoforma. The unmodified-CDS baseline accepts 41,592 exact translations, flags 17,399 non-triplet lengths and records 125 translation mismatches. These immutable strict results are retained separately.

The genome/GFF audit reconstructs transcription-ordered CDS blocks, checks assembly bounds, strand, non-overlap and every internal phase. It accepts two explicitly recorded published export conventions: the full annotated genomic CDS span, or that exact span with its annotated initial phase bases already removed. It derives complete codons from the genomic span, applying the initial phase once and omitting only a terminal 1–2-base remainder. It never searches reading frames or changes genetic codes to obtain a match. Omitted bases and original genomic segments are recorded; terminal stop codons remain flagged in the DNA output.

The first projection implementation required the full genomic span to equal the published string and flagged 12,431 Pirum/Abeoforma records. Every one is explained by the second export convention: 7,097 Pirum and 5,334 Abeoforma sequences already omit the annotated initial phase bases. This is a source representation difference, not a demonstrated genome sequence error. Version 1 results remain archived; producer source is preserved in commit 0148d61.

Version 2 verifies **all 59,116 genome-linked complete-codon spans**: 7,535 Corallochytrium, 12,463 Chromosphaera, 21,835 Pirum and 17,283 Abeoforma. Partial boundaries occur in 25, 1,020, 13,565 and 9,757 records respectively. The 125 strict Chromosphaera translation mismatches are explained by annotated nonzero initial phases. All translated protein residues must match after removing at most one terminal stop from the translation.

```bash
python scripts/validate_published_outgroup_cds.py --output results/cds/published-outgroup-translation-v1
python scripts/project_published_outgroup_codons.py --strict results/cds/published-outgroup-translation-v1 --output results/cds/published-outgroup-codon-projection-v2
python -m unittest discover -s tests -p test_published_cds_translation.py
python -m unittest discover -s tests -p test_codon_projection.py
```

Use new output directories for reruns. Seven focused tests cover strict stop/partial-codon handling and genomic projection, including reverse strand and avoiding double phase removal. Independent readback checks all artifact hashes and reproduces all 59,116 translations; the 12,431 version-1 exceptions exactly equal the version-2 pretrimmed-export set. Receipts and readback evidence are versioned under `metadata/published_outgroup_*`; large per-protein audits and nucleotide FASTAs remain in the corresponding `results/cds` directories. The historical external source inventory records the acquisition state before these checks; these newer receipts supply downstream validation.

Complete-codon spans do not establish complete genes or suitability for selection inference. Partial-boundary inclusion policies, gene-representative links, codon alignments, genetic-code review for additional taxa and family-specific divergence assessment remain necessary.

## Full NCBI strict translation audit running

`audit_ncbi_cds_translation.py` now audits all 519 NCBI-backed taxa against their exact assembly-matched GFF and normalized protein input. Every source hash is checked. It associates CDSs through exact versioned protein IDs, reports missing IDs and proteins, and excludes multiply linked CDS records from direct-match classification. GFF `transl_table` values supply the translation code; conflicts remain exceptions. Missing explicit codes use a clearly recorded table-1 assumption for this diagnostic, requiring subsequent code review before selection eligibility.

The audit translates unmodified triplet-length DNA, removes at most one terminal stop from the translated string, and requires every protein residue to match. Non-triplet lengths, different translations and GFF exception/recoding/pseudogene flags remain explicit. It does not trim phases, repair initiation residues, implement translational exceptions or search codes. Partial annotation flags remain attached to direct matches, so an exact translation is not evidence of a complete gene. Raw CDS headers, locations and annotation flags are retained in per-taxon audit tables.

```bash
python scripts/audit_ncbi_cds_translation.py --output results/cds/ncbi-strict-translation-v1
python -m unittest discover -s tests -p test_ncbi_cds_translation.py
```

Execution uses one worker on the existing host, with a 0.2–4-hour scheduling estimate, 4 GB memory and 5 GB output allowance recorded before launch in `metadata/ncbi_cds_translation_resource_plan.json`. Per-taxon completed receipts support restart only after source, configuration and output hash verification; incomplete taxa are recomputed. A single-writer lock prevents simultaneous producers. Four tests cover alternative CUG translation, missing/conflicting code provenance, no partial-codon or initiator repair, and at-most-one terminal stop removal. The run log is `logs/ncbi_strict_cds_translation_v1.log`; the full receipt is written only after all taxa finish. This stage is running, not completed; codon alignments, representative-gene integration and selection tests remain pending.

## Completed full marker CDS identity index

The full marker index links all 59,840 existing marker/taxon records to available CDS sources by exact protein identifier. NCBI CDS FASTAs are read unchanged; the four published outgroups use the verified genome-projected codon spans, while the other three external taxa use their verified extracted subsets. Every source CDS, normalized protein and representative-decision file is hash checked, and each marker protein sequence is rechecked against its original sequence hash.

Duplicate CDS records for a protein remain ambiguous even if their DNA strings are identical. Missing CDSs remain explicit. Unique records are exported under marker/taxon identifiers, accompanied by source record IDs, nucleotide hashes, source conventions, gene IDs and representative decisions. Alternative products are retained as flagged marker identities rather than replaced with another isoform. This is a source index, not translation qualification: NCBI records must be joined to the completed strict audit, and external partial-boundary policies must accompany codon-alignment inclusion.

```bash
python scripts/index_marker_cds_sources.py --output results/cds/marker-source-index-v1
python -m unittest discover -s tests -p test_marker_cds_index.py
```

The immutable index completed with one worker; a 2–30-minute scheduling estimate, 2 GB memory and 1 GB output allowance were recorded in `metadata/marker_cds_index_resource_plan.json`. Two tests cover missing/unique records and duplicate ambiguity. Execution is logged in `logs/marker_cds_source_index_v1.log`. Large FASTA and per-marker tables remain outside Git.

All **59,840 marker/taxon links across 526 taxa** have exactly one associated source CDS. The exported FASTA contains 100,455,973 nucleotides. Independent readback checks the complete marker identity universe, every original marker provenance field, source-record multiplicity, all exported sequence hashes/lengths/alphabet and the exact 526-taxon manifest membership.

Gene decisions remain important despite complete CDS availability: 59,523 markers are selected representatives of unique genes, four are alternative products, and 313 retain unresolved/provisional gene mappings (124 unmapped and 189 provisional ORFs). These are not silently removed or substituted. Source indexing does not yet qualify all NCBI translations or establish complete genes.

```bash
python scripts/verify_marker_cds_index.py --index results/cds/marker-source-index-v1 --output results/cds/marker-source-readback-v1
```

Receipts, full taxon coverage and the 317 marker gene/isoform exception rows are versioned as `metadata/marker_cds_source_index_receipt.json`, `metadata/marker_cds_source_readback_receipt.json`, `metadata/marker_cds_taxon_coverage.tsv` and `metadata/marker_cds_gene_exceptions.tsv`. The next codon-alignment gate joins exact CDS identities to completed translation results and external boundary flags, retaining genetic-code provenance.

## NCBI code defaults and marker boundary audit

NCBI documents table 1 as the default when a translation-table qualifier is absent; nonstandard values are supplied explicitly. Its GFF documentation also describes source-region translation codes and strand-aware partial CDS boundaries. These conventions provide source-format evidence rather than choosing a code because its translation matches. See [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/) and [NCBI GFF3 format](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/) (accessed 2026-09-13).

The already running full-proteome audit retains its original `table_1_assumption_no_explicit_code` labels and immutable producer. The new marker-specific boundary audit records the documented default as distinct from explicit CDS and source-region codes. It checks disagreements and retains the original source attributes and ordered CDS segments. It separately records initial phase, internal phase consistency, overlaps, partial 5-prime and 3-prime ends, internal partial boundaries, pseudogene/translation-exception flags and strict unmodified-CDS translation. Direct translation matches do not override annotation issues.

```bash
python scripts/audit_ncbi_marker_boundaries.py --output results/cds/ncbi-marker-boundaries-v1
python -m unittest discover -s tests -p test_marker_cds_boundaries.py
```

This full 519-taxon audit is running with one worker and source-hash checks. Prelaunch scheduling estimate: 3–30 minutes, 2 GB memory and 1 GB output; see `metadata/marker_boundary_resource_plan.json`. Five tests cover reverse-strand partial starts, internal versus terminal boundaries, regional/default codes, code conflicts and split-codon phase continuity. The output supplements the marker source index and full-proteome translation audit. It does not reconstruct NCBI genomic sequence, repair frames or establish selection eligibility. External marker boundaries retain their separate verified extraction/projection policies. Producer log: `logs/ncbi_marker_boundaries_v1.log`.
