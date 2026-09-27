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

## Completed full NCBI strict translation audit

`audit_ncbi_cds_translation.py` now audits all 519 NCBI-backed taxa against their exact assembly-matched GFF and normalized protein input. Every source hash is checked. It associates CDSs through exact versioned protein IDs, reports missing IDs and proteins, and excludes multiply linked CDS records from direct-match classification. GFF `transl_table` values supply the translation code; conflicts remain exceptions. Missing explicit codes use a clearly recorded table-1 assumption for this diagnostic, requiring subsequent code review before selection eligibility.

The audit translates unmodified triplet-length DNA, removes at most one terminal stop from the translated string, and requires every protein residue to match. Non-triplet lengths, different translations and GFF exception/recoding/pseudogene flags remain explicit. It does not trim phases, repair initiation residues, implement translational exceptions or search codes. Partial annotation flags remain attached to direct matches, so an exact translation is not evidence of a complete gene. Raw CDS headers, locations and annotation flags are retained in per-taxon audit tables.

```bash
python scripts/audit_ncbi_cds_translation.py --output results/cds/ncbi-strict-translation-v1
python -m unittest discover -s tests -p test_ncbi_cds_translation.py
```

Execution uses one worker on the existing host, with a 0.2–4-hour scheduling estimate, 4 GB memory and 5 GB output allowance recorded before launch in `metadata/ncbi_cds_translation_resource_plan.json`. Per-taxon completed receipts support restart only after source, configuration and output hash verification; incomplete taxa are recomputed. A single-writer lock prevents simultaneous producers. Four tests cover alternative CUG translation, missing/conflicting code provenance, no partial-codon or initiator repair, and at-most-one terminal stop removal. The run log is `logs/ncbi_strict_cds_translation_v1.log`; the full receipt is written only after all taxa finish. This stage is complete. Marker codon projection and independent readback subsequently completed as described below; full-family codon preparation and selection tests remain pending.

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

This full 519-taxon audit completed with one worker and source-hash checks. Prelaunch scheduling estimate: 3–30 minutes, 2 GB memory and 1 GB output; see `metadata/marker_boundary_resource_plan.json`. Five tests cover reverse-strand partial starts, internal versus terminal boundaries, regional/default codes, code conflicts and split-codon phase continuity. The output supplements the marker source index and full-proteome translation audit. It does not reconstruct NCBI genomic sequence, repair frames or establish selection eligibility. External marker boundaries retain their separate verified extraction/projection policies. Producer log: `logs/ncbi_marker_boundaries_v1.log`.

## Codon projection onto fixed protein alignments

The codon projection implementation preserves the 125 full MAFFT marker alignments and their previously hashed 50-percent protein occupancy masks. It verifies full-protein translation before selecting columns, inserts complete triplet gaps for protein gaps, and removes at most one verified terminal stop codon. A mismatch outside the retained columns still fails projection. Output columns map explicitly to the original protein alignment.

NCBI sequences require strict translation agreement and consistent annotation metadata; nonzero initial phases, code conflicts and translational exceptions remain excluded pending coordinate review. Exact partial translations can enter these diagnostic alignments with boundary flags. External sequences inherit their verified extraction/projection policies. Gene alternatives and unresolved mappings remain flagged. Genetic codes are recorded per sequence; mixed-code marker alignments cannot be passed indiscriminately to a single-code selection model.

```bash
python scripts/build_marker_codon_alignments.py --output results/cds/marker-codon-alignments-v1
python scripts/verify_marker_codon_alignments.py --alignments results/cds/marker-codon-alignments-v1 --output results/cds/marker-codon-readback-v1
python -m unittest discover -s tests -p test_codon_alignment_projection.py
```

Five projection tests cover masked/gapped correspondence, nonstandard codes, mismatches outside retained columns, stop/partial-codon rejection and invalid masks. The independent verifier reads every output codon back against original masked protein columns and checks the complete identity/status universe. The prelaunch plan allows one worker, 3 GB memory, 1 GB output and 1–20 minutes for projection. These are diagnostic codon alignments; clade-specific saturation, copy/gene policies, alignment sensitivity and selection-model eligibility remain separate requirements. Execution/completion status is recorded in the progress log and result receipts.

The completed boundary audit covers **59,269 NCBI marker sequences**: 58,967 strict exact translations, 280 non-triplet CDS lengths and 22 translation mismatches. All have consistent ordered CDS feature metadata under the checked strand/phase/overlap rules. Partial 5-prime and 3-prime boundaries occur in 919 and 965 records respectively (overlap is possible); no internal partial boundary was found. There are 201 nonzero initial phases. These flags remain independent of strict translation status.

The code inventory comprises 58,562 table-1 defaults and 707 explicit CDS codes: 586 table-12, 120 table-26 and one table-6 marker sequences. No marker required source-region code fallback or showed a code conflict. Thus documented defaults resolve code provenance for this marker set without fitting codes to amino-acid agreement. Metadata consistency does not prove that all gene models are biologically correct.

`metadata/ncbi_marker_boundary_receipt.json` preserves source and output provenance; `metadata/ncbi_marker_boundary_readback.json` records independent artifact/identity checks and recounted totals. A compact 1,673-row review table retains translation, boundary or annotation flags in `metadata/ncbi_marker_boundary_review.tsv`; full ordered genomic segments remain outside Git. Codon projection has been launched after successful boundary completion; its final receipt and independent codon readback determine completion.

Projection completed for **all 125 markers**, accepting **59,334 marker/taxon sequences** and excluding 506 with explicit reasons. The unchanged masks retain 63,750 codon columns across markers. Each marker has 418–502 accepted taxa; 123 markers contain more than one translation code. The exclusions comprise 267 non-triplet sequences, 188 initial-phase review cases, 16 annotation-exception cases, 13 with both phase and non-triplet flags, 21 translation mismatches and one annotation-exception/mismatch case. These counts are mutually exclusive status combinations.

The producer receipt is `metadata/marker_codon_alignment_receipt.json`. Compact exclusion/boundary/gene review rows are in `metadata/marker_codon_alignment_review.tsv`. The separate verifier completed successfully: every one of 28,328,533 non-gap codons in 59,334 sequence rows translated under its recorded code to the exact original masked amino-acid column. All artifact hashes, identities, gap triplets and marker summary totals passed. Evidence is in `metadata/marker_codon_readback_receipt.json` and `metadata/marker_codon_summary.tsv`. Log: `logs/marker_codon_readback_v1.log`. No selection model has been fitted to these codon alignments.

## Coverage screen for groups sharing genus labels and codes

`assess_codon_group_coverage.py` requires the completed independent codon readback. It evaluates all 125 markers for fungal genus labels represented by at least four manifest entries, splitting each marker/group by translation code. These labels define candidate groups for assessment, not inferred clades. Entries are not assumed to be unique species or independent ecological transitions.

Two policies retain either all translation-aligned records or only records without the existing annotation, gene-representative and taxon-label flags. Both use the unchanged protein-derived column mask. Complete unambiguous codons define called data. The diagnostic coverage screen asks whether at least four entries remain with at least 100 columns called in at least 80 percent of that subset; it also reports fully called and variable covered columns. Empty groups and excluded entries remain explicit.

```bash
python scripts/assess_codon_group_coverage.py --output results/cds/genus-code-coverage-v1
python -m unittest discover -s tests -p test_codon_group_coverage.py
```

Three tests verify ambiguity/gap treatment, variable covered columns and empty subsets. The resource plan allows one worker, 1 GB memory, 0.1 GB output and 1–10 minutes. Passing this coverage screen does not establish selection eligibility: supported group phylogenies, gene/copy reconciliation, alignment sensitivity, divergence/saturation, taxon identity and model adequacy still require assessment. This screen does not estimate dN/dS or test ecological effects.

The completed screen covers **16 fungal genus labels and 113 manifest entries**, producing 4,168 marker/code/policy rows. There are 1,826 passing marker/code groups under all-aligned inclusion and 1,712 after excluding recorded flags. All 16 labels have passing groups under inclusive screening; 15 retain some under the stricter policy. The stricter passing-marker counts range from zero (Serendipita) to 125.

Serendipita illustrates a taxon-identity dependency: its five entries include three unnamed species records already flagged by the label audit. Excluding those leaves only two named entries, below the four-entry threshold; the change is not evidence of poor sequence coverage or biological absence. Candida contributes multiple code groups and remains a label-defined group requiring phylogenetic review. Coverage passing is not permission to treat these taxa as a resolved clade or replicated ecological transitions.

All output hashes, unique group identities, retained-taxon counts, stricter-policy subset membership and threshold totals passed readback checks. This readback does not independently recompute every coverage base. Summaries, exact genus membership and receipts are in `metadata/codon_group_coverage_*` and `metadata/codon_group_membership.tsv`; the full row table remains in `results/cds/genus-code-coverage-v1`.

## Completed full-proteome translation results and reconciliation

The complete audit covers **519 NCBI taxa, 5,847,336 CDS records and 5,843,347 normalized proteins**. Every normalized protein has a source CDS association, although an association does not guarantee a unique or valid translation. The mutually exclusive CDS classifications are:

| Classification | CDS records |
| --- | ---: |
| Exact translation | 5,752,457 |
| Non-triplet length after earlier gates | 85,683 |
| Translation mismatch | 2,491 |
| Annotation exception requiring review | 2,715 |
| Missing/ambiguous protein identifier | 3,988 |
| Multiple CDS records for one protein | 2 |

Classification order matters: counts above are status partitions, not independent counts of every possible overlapping flag. The acquisition audit's raw non-triplet count therefore need not equal the classified non-triplet total. All DNA and raw headers remain preserved; no frame shifts, initiator repairs or exception recoding were introduced.

The full audit records 5,793,946 default-table-1 evaluations under its original assumption label, and explicit tables 6 (13,180 records), 26 (5,666), 12 (29,979), 3 (57), 4 (493) and 5 (25). Rows rejected before code assignment are excluded from these counts. The marker-specific review verifies documented default and CDS/region-code provenance for the marker subset; source-region fallback was not implemented in this historical full-proteome producer and remains a consideration for later nonmarker qualification.

```bash
python scripts/summarize_ncbi_cds_audit.py --audit results/cds/ncbi-strict-translation-v1 --output results/cds/ncbi-strict-readback-v1
```

The completed readback verifies source-inventory, configuration and producer hashes, all per-taxon output hashes, unique CDS identities, full 519-taxon coverage and every status/code total. It reconciles all 59,269 NCBI marker protein codes with the independently executed boundary audit: 59,252 have identical strict translation statuses, while 17 differences arise because the full audit stops at an annotation-exception gate before translating. There are no unexplained disagreements. The compact 17-row table is `metadata/ncbi_marker_translation_gate_differences.tsv`.

Full receipts, the readback and taxon summaries are versioned as `metadata/ncbi_cds_translation_receipt.json`, `metadata/ncbi_cds_translation_readback.json` and `metadata/ncbi_cds_taxon_translation_summary.tsv`. Per-CDS tables remain in `results/cds/ncbi-strict-translation-v1`. This readback does not retranslate every nonmarker sequence independently or validate its genomic coordinates. Exact translation still does not establish selection eligibility, gene-model correctness or complete biological genes.

## Observed divergence within coverage-screened groups

The completed observed-difference screen evaluates 3,538 marker/genus/code/policy groups, yielding 84,796 pair/policy rows and **43,235 distinct marker/code/taxon pairs**. There are 43,194 pairs with at least 100 shared called codons under inclusive selection and 41,529 under the stricter recorded-flag policy. Low-overlap rows remain explicit. The same pair repeated under two policies has identical metrics and is not a replicate.

Every metric uses pairwise shared, unambiguous complete codons. Outputs count different codons, different amino acids, different codons encoding the same amino acid, and differences at each nucleotide position. Fractions are uncorrected observations. Same-amino-acid codon differences are not reconstructed synonymous substitutions; no dS, dN/dS, saturation test or branch rate is estimated.

```bash
python scripts/assess_codon_pair_divergence.py --output results/cds/genus-code-pair-divergence-v1
python scripts/plot_codon_group_divergence.py --input results/cds/genus-code-pair-divergence-v1 --output results/cds/genus-code-divergence-figure-v1
python -m unittest discover -s tests -p test_codon_pair_divergence.py
```

Four tests cover code-dependent translation, shared-site denominators, synonymous observed codon differences, zero shared data and stop rejection. Readback checks source/output hashes, count partitions and bounds, fraction denominators, overlap thresholds and identical repeated metrics. It does not independently recompute all pairs from DNA. Source receipts and group summaries are versioned under `metadata/codon_pair_divergence_*` and `metadata/codon_group_pair_summary.tsv`; large pair tables remain outside Git.

![Observed codon divergence across marker groups](figures/codon_group_divergence.svg)

The figure shows 1,712 stricter-policy marker/group medians across 15 fungal genus labels. It is descriptive: pairs/points share taxa, sites and ancestry, and genus labels are not verified clades. The SVG and PNG were rendered and visually inspected; figure provenance is `metadata/codon_group_divergence_figure_receipt.json`.

The 20 largest marker/group median amino-acid differences are retained for **alignment/orthology review**, not designated accelerated evolution (`metadata/codon_divergence_alignment_review.tsv`). The largest is marker 4986044at2759 in Aspergillus: ten entries, at least 232 shared called codons per pair, median amino-acid difference 0.7371 and median third-position difference 0.6767. Group-specific realignment, existing profile-alignment sensitivity and protein/domain correspondence need review before biological interpretation. Other leading candidates include marker 776280at2759 in Naganishia and Tilletia. Coverage and translation alone do not resolve these concerns.

## Targeted alignment and domain review

All 20 leading divergence cases were realigned within their recorded genus/code subsets using the exact original full protein sequences and MAFFT auto with one thread. Residue identity preservation passed for every output. The original masked pairwise divergences were independently reproduced from amino-acid alignments while restricting comparisons to residues with canonical source codons. Three focused tests cover residue versus column coordinates, source-codon ambiguity and mask offsets.

```bash
python scripts/realign_codon_review_groups.py --output results/cds/outlier-realignment-v1
python scripts/audit_aspergillus_tfiib_marker.py --output results/cds/aspergillus-domain-hit-review-v1
python -m unittest discover -s tests -p test_residue_pair_correspondence.py
```

The control covers 20 selected cases and 163 taxon pairs; it is targeted diagnosis of observed outliers, not a representative calibration or separate pilot. Within-group MAFFT can choose a different algorithm from full-dataset MAFFT, and its 80-percent occupancy mask differs from the original 50-percent global mask. Changes therefore combine alignment, taxon context and selected sites. Across these cases, the local-minus-global median amino-acid difference ranges from −0.0371 to +0.0522. None of the high differences simply disappears. Exact residue-pair retention and Jaccard statistics make correspondence sensitivity explicit.

For Aspergillus marker **4986044at2759**, the group median changes from 0.7371 to 0.7605; median retention of original residue pairs is 0.6932 and median correspondence Jaccard is 0.4926. Hash-verified Pfam evidence separates six 713–752-residue proteins with BRF1-family hits from four 350–351-residue proteins without those hits. Both subsets have TFIIB hits; the shorter subset also has Zn_Ribbon_TF hits. These are observed HMM annotations, not validated functions or resolved architectures.

The 45 pairs divide into 24 between-subset comparisons, 15 within the BRF1-hit subset and six within the other subset. Global masked median amino-acid differences are **0.7500 between subsets**, **0.1336 within BRF1-hit proteins**, and **0.0460 within proteins without BRF1 hits**. Local-alignment values are 0.7645, 0.0581 and 0.0808 respectively. This exploratory partition explains why between-type contrasts dominate the pooled median, while leaving the origin of those types unresolved.

Mixed homologous protein types or hidden paralogy is now a concrete competing explanation. Gene-family phylogeny, broader homolog sampling and reconciliation are required before this group can support within-ortholog acceleration or selection claims. Original BUSCO markers and guide-tree inputs remain preserved for sensitivity comparisons; no domain gain/loss or validated function is inferred. Case artifacts, sequence/identifier correspondence and all realignment outputs passed hash readback. Versioned receipts/summaries are `metadata/codon_outlier_realignment_*`, `metadata/aspergillus_marker_*` and `metadata/codon_outlier_case_readback.json`; detailed alignments, pair correspondences and HMM hit records retain source provenance.

## Broad homolog recovery for the mixed-domain case

The protected completed OrthoFinder core table was checked against the checksum recorded before full assignment. Only one of the ten focal Aspergillus taxa is in that 64-taxon core: F746128, whose selected protein XP_754348.1 belongs to OG0001567 (62 core taxa). The other nine are outside the core, not demonstrated homolog absences. This reference association cannot resolve the observed contrast. Exact membership and source hashes are in `metadata/aspergillus_marker_core_trace.json`.

A focused full-representative search has completed for Pfam profiles PF00382.25 (TFIIB), PF07741.19 (BRF1) and PF08271.18 (Zn_Ribbon_TF). It searches all 5,654,720 additional unique proteins from the 526-taxon design; existing marker hits are available for subsequent combination using the full representative-protein mappings. This recovers alternative homolog candidates beyond the single BUSCO-selected protein without waiting for every full-Pfam profile to finish.

```bash
python scripts/search_tfiib_family_domains.py --output results/domains/tfiib-family-search-v1
```

The controller verifies the complete Pfam HMM library, all full search-input artifacts, reusable marker hits and source compatibility; extracts the exact three accession-version profiles; and runs the same pinned HMMER 3.4 binary with sequence/domain gathering thresholds, two workers within one search process, a coordinator and seed 42. It records profile, binary, input and command hashes and requires a completed report containing all three queries before issuing a final receipt. No full-family result is claimed while the search is active. The prelaunch plan allows 2–120 minutes, 8 GB memory and 1 GB output on the existing host. Profiles/configuration are under `results/domains/tfiib-family-search-v1`; versioned config and plan are `metadata/tfiib_family_search_config.json` and `metadata/tfiib_family_search_resource_plan.json`.

HMM hits nominate homolog/domain candidates; they do not resolve orthogroups, validate functions, demonstrate loss or reconstruct duplication history. Combining sequence hits with all gene-representative identities, assessing domain correspondence and constructing/reconciling family trees remain necessary. E-values from marker and additional-protein searches have different target database sizes; common gathering thresholds are the inclusion criterion.

The focused search completed in **64.54 seconds**, returning **2,109 domain-hit rows** among the 5,654,720 additional sequences. All three query profiles completed, and output/producer hashes and row/profile counts passed readback. `metadata/tfiib_family_search_receipt.json` and `metadata/tfiib_family_search_readback.json` record the measured result. These are additional-sequence hits; combination with existing marker hits and mapping back to all representative genes remain the next steps before family alignment or reconciliation.

## Full candidate identity collection completed

Combined all 2,109 focused additional-protein hits with 1,649 reusable marker hits from the same three Pfam profiles. Mapping through the full representative-protein input restores **1,485 protein/gene entries (1,434 unique sequences) across 524 of 526 taxa**. Every hit-bearing sequence has a representative association. Candidate gene annotations comprise 1,477 unique-gene mappings, six provisional ORFs and two unmapped entries. No sequence-similarity deduplication removes gene/taxon identity from the exported family inputs.

```bash
python scripts/collect_tfiib_family_candidates.py --output results/domains/tfiib-candidates-v1
```

The collector verifies both search partitions and their compatible input/annotation provenance, streams full sequence/protein mappings, checks every selected sequence hash and domain coordinate bound, and restores gene-representative decisions. A separate readback verified every exported sequence/hash/length and all unique protein identities. All source hits remain available, with their search partition explicit.

Scedosporium apiospermum (F563466) and Abeoforma whisleri (OFS5426458) have no passing candidate in this three-profile screen. This is not evidence of biological absence or gene loss; divergence, annotation quality and profile sensitivity remain alternatives. Candidate proteins are a union of domain matches, not one established homologous full-length family. There are 1,099 proteins with two TFIIB hits, 340 with one, 44 with none and two with three/four hits. Domain-copy correspondence and coverage must be considered before alignment.

All ten focal Aspergillus genomes contain **two candidates on distinct uniquely annotated genes**, one with a BRF1 hit and one without. In each genome exactly one of those candidates is the original marker 4986044at2759. The marker chooses the BRF1-hit type in six taxa and the other type in four, while both types are recovered in every focal genome. This makes inconsistent copy sampling a concrete competing explanation for the pooled marker divergence. It does not establish duplication timing, ortholog relationships or a domain gain/loss event; those require family trees and reconciliation.

Versioned outputs include `metadata/tfiib_candidates_receipt.json`, `metadata/tfiib_candidate_protein_mapping.tsv`, `metadata/tfiib_candidate_taxon_coverage.tsv`, `metadata/tfiib_candidates_readback.json` and the 20 focal copies in `metadata/aspergillus_family_candidate_copies.tsv`. Full candidate FASTA and domain coordinates remain under `results/domains/tfiib-candidates-v1`.

## Ordered domain-pair alignment and supported tree run

The completed candidate screen selected proteins with exactly two nonoverlapping TFIIB GA hits, each covering at least 70 percent of PF00382.25. Domain count, overlap and coverage exceptions remain explicit; it does not choose the best two hits from a protein with additional copies. Of 1,485 candidates, **1,040 proteins from 512 taxa** enter the paired alignment; 386 have a different hit count, 58 fail the per-domain coverage threshold and one has overlapping hits.

```bash
python scripts/align_tfiib_repeat_pairs.py --output results/phylogeny/tfiib-domain-pairs-v1
python -m unittest discover -s tests -p test_tfiib_pair_selection.py
python scripts/run_tfiib_family_tree.py --output results/phylogeny/tfiib-domain-tree-v1
```

The two domain segments are ordered N-to-C, separately aligned against the same pinned HMM with HMMER 3.4, and concatenated as two 92-match-state blocks (**184 columns**). Full Stockholm alignments preserve every input domain residue. Every retained match-state residue maps explicitly to its original full-protein coordinate. Independent readback verified all 191,360 mapping cells and all 175,699 non-gap residues. Three focused tests cover positional ordering, refusal to choose two of three hits, overlap and insufficient coverage.

This positional repeat correspondence is a working hypothesis; repeat-specific phylogenies, domain order and possible gene conversion require later sensitivity checks. It does not establish full-length orthology or a rooted duplication history. All excluded proteins remain in the candidate inventory.

A supported unrooted gene-tree search is now running on the full 1,040-protein paired-domain alignment. IQ-TREE 3 uses restricted model selection among LG/WAG/JTT with empirical frequencies and gamma rates, 1,000 SH-aLRT replicates and 1,000 ultrafast-bootstrap replicates with NNI refinement; bootstrap trees are retained. Four threads and 8 GB memory are allocated with a prelaunch 0.25–8-hour scheduling allowance. Configuration, binary/input hashes and seed 20260913 are recorded. The controller accepts completed output only after checking all original tip identities, finite nonnegative branches and 1,000 bootstrap trees. Identical inputs or zero-length branches must not be interpreted as resolved duplication events. Checkpoint resume requires unchanged configuration and source hashes.

Versioned evidence: `metadata/tfiib_domain_alignment_receipt.json`, `metadata/tfiib_domain_alignment_readback.json`, `metadata/tfiib_domain_pair_candidate_audit.tsv`, `metadata/tfiib_family_tree_run_config.json` and resource plans. Tree execution log: `results/phylogeny/tfiib-domain-tree-v1/stdout.log`. Tree completion, branch-support interpretation and species-tree reconciliation are not yet claimed.

## Repeat-specific phylogenetic sensitivity

The verified paired-domain matrix was split into its original two 92-column blocks with the same 1,040 gene identities and no additional filtering. The first repeat has 862 distinct aligned sequences, 91 variable/parsimony-informative columns and one constant column. The second has 894 distinct sequences, 90 variable/parsimony-informative columns, one constant and one all-gap/ambiguous column. These descriptive counts do not establish sufficient phylogenetic information or independence.

```bash
python scripts/prepare_tfiib_repeat_sensitivity.py --output results/phylogeny/tfiib-repeat-sensitivity-inputs-v1
python scripts/run_tfiib_repeat_trees.py --output results/phylogeny/tfiib-repeat-trees-v1
```

Independent readback confirms that all repeat identities match the original matrix and that recombining every pair exactly reproduces the original 184-column sequence. Input receipt, per-gene coverage and readback are versioned under `metadata/tfiib_repeat_input_*` and `metadata/tfiib_repeat_coverage.tsv`.

Both repeat-specific IQ-TREE searches are running, each with the same restricted LG/WAG/JTT + empirical-frequency/gamma model set, 1,000 SH-aLRT replicates and 1,000 NNI-refined ultrafast-bootstrap replicates. They use two threads and a 4 GB memory limit per job, with distinct recorded seeds. The combined-domain search continues separately. Source/configuration hashes and exact input identities are checked; the controllers require valid completed trees and all 1,000 bootstrap outputs before claiming success. Metadata records both run configurations and the resource allowance.

The comparison holds gene sampling fixed while changing which repeat supplies the sites. Discordance could reflect limited information, model/alignment issues, repeat history or other biological processes; it will not by itself prove gene conversion. Domain-pair and repeat-specific tree completion, support-aware concordance and reconciliation remain pending.

## Genus-label checks against completed gene trees

A frozen, explicitly incomplete snapshot contains 17 of 125 marker trees,
with 8,031 audited internal edges. `assess_genus_tree_splits.py` evaluates the
16 genus labels used in the codon coverage screen against this exact snapshot.
For each marker it records missing genus members and requires at least four
observed members and two other taxa. It checks whether the observed group is
one side of an unrooted internal edge. An incompatible edge must have all four
intersections of the two bipartitions nonempty; the result is independent of
which side is represented first.

Of 272 marker/genus rows, 250 are assessable and 20 contain an incompatible
edge with SH-aLRT support at least 80. Naganishia has such conflicts in five of
ten assessable markers (129234at2759, 260326at2759, 340246at2759,
345792at2759 and 4747214at2759). The remaining five have a separating edge.
These are descriptions of gene-tree splits, not a rooted species-monophyly test,
independent ecological transitions, orthology confirmation or a selection
result. SH-aLRT is not a posterior/ bootstrap probability; early-finishing genes
may not represent the full set. Alternative causes include paralogy, alignment
or inference error and biological discordance. No genus or marker is removed
from the source data by this screen.

```bash
python scripts/audit_marker_tree_support.py \
  --trees results/phylogeny/marker-gene-trees-v2 \
  --output results/phylogeny/marker-tree-support-v3 --allow-incomplete
python scripts/assess_genus_tree_splits.py \
  --audit results/phylogeny/marker-tree-support-v3 \
  --trees results/phylogeny/marker-gene-trees-v2 \
  --output results/phylogeny/genus-tree-splits-v1
```

Use new outputs for reruns. Two split-compatibility tests pass; output hashes
and the full 17-by-16 row universe were checked. Per-marker/group rows and
summary/source receipts are versioned. Review these conflicts with the full
marker set, family assignments and the species framework before interpreting
genus-based codon-model results.


## Within-genus diagnostic inputs

Prepared all 2,084 previously strict genus/code groups with current curated-hybrid
exclusions. The full ledger retains 372 insufficient-coverage cases and 1,712
groups ready for tree/divergence diagnostics. Ready groups have 4–10 taxa and
100–2,342 codons, with 5,728,643 observed taxon–codon cells. Of these groups,
1,613 use NCBI code 1 and 99 use code 12. Fourteen ready groups retain the TFIIB
copy-reconciliation caveat; none is declared suitable for selection testing.

Within each group, retain source codon columns called in at least 80% of retained
taxa, including invariant columns. Missing/noncanonical codons become whole
`???` triplets with matching `?` protein states. Export exact mappings to the
original codon and MAFFT protein columns. This preserves existing homology
assumptions; it does not independently validate the alignment. Current hybrid
exclusions affected membership in 14 cases without changing the total passing
coverage count.

`scripts/prepare_genus_codon_diagnostics.py` creates immutable paired DNA/protein
inputs in `results/cds/genus-codon-diagnostic-inputs-v1`.
`scripts/readback_genus_codon_diagnostics.py` independently passed every case,
updated membership, occupancy mask, exported codon, translation and source
column. Receipts and full case summary are versioned under
`metadata/genus_codon_diagnostic_*`.

Next stages are group-specific phylogenies and codon-model divergence diagnostics,
with alignment, gene-copy, recombination and genetic-code review. A genus label
does not establish monophyly; tree availability or a coverage pass does not prove
absence of synonymous saturation or validate a selection claim.


## Supported within-genus nucleotide trees

The reproducible information screen retains 1,655 of the 1,712 prepared groups
for supported tree searches: at least four distinct aligned nucleotide strings
and one parsimony-informative nucleotide column (at least two canonical states
observed twice each). Missing states are ignored when counting informative
columns, but retained in aligned strings. This is a workflow screen, not proof
of phylogenetic adequacy. The other 57 cases remain recorded for possible
analyses with externally justified topologies.

All 1,655 searches are running with IQ-TREE 3.0.1, unpartitioned GTR+F+G4 DNA,
1,000 SH-aLRT and 1,000 ultrafast bootstrap replicates with bootstrap NNI;
bootstrap trees are retained. Four concurrent workers use one CPU and a 2 GB
memory allowance each. Resource planning reserves 10 GB output and 1–96 hours
on the existing host; this range is not a measured completion bound. Case seeds,
alignments, executable, scripts and configuration are checksum-pinned. The
runner preserves genetic-code and gene-copy caveats for later codon models.
Execution checks exact tip sets and finite nonnegative edges; independent full
report/model/support audit and codon divergence diagnostics remain pending.

The information table was regenerated and matched byte for byte across all
1,712 cases. At the timestamp recorded in
`metadata/genus_codon_tree_execution_observation.json`, 287 cases had execution
receipts; this is a running observation, not full completion.

```bash
python scripts/screen_genus_codon_tree_information.py \
  --inputs results/cds/genus-codon-diagnostic-inputs-v1 \
  --output /tmp/genus_codon_tree_information.tsv
# Compare the regenerated screen to metadata/genus_codon_tree_information.tsv.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/run_genus_codon_trees.py \
  --inputs results/cds/genus-codon-diagnostic-inputs-v1 \
  --readback metadata/genus_codon_diagnostic_readback.json \
  --information metadata/genus_codon_tree_information.tsv \
  --resources metadata/genus_codon_tree_resource_plan.json \
  --output results/phylogeny/genus-codon-trees-v1
```

Do not start a second instance while this queue is active. Completed cases can
be reused only with matching configuration and artifact hashes. These unrooted,
within-genus/code trees do not establish genus monophyly, orthology, species-tree
relationships, absence of synonymous saturation, or selection eligibility.


`scripts/audit_genus_codon_trees.py` independently parses every completed ML
and bootstrap tree with Bio.Phylo. It checks the exact taxon grid, finite
nonnegative branches, resolved unrooted splits, input/configuration/artifact
hashes, recorded model/support settings, finite likelihood and gamma shape,
and the reported total tree length. It recounts each ML split across all 1,000
saved bootstrap trees; reported integer UFB support must agree within rounding
(0.5 percentage points). SH-aLRT ranges are checked but its statistics are not
independently recomputed. Warnings and near-zero/very long edges are retained.

Use `--allow-incomplete` only for a new immutable snapshot of completed cases;
the receipt lists all pending cases. Without that option, the full batch receipt
and exact case-receipt grid are required. Neither mode proves model adequacy,
optimization convergence or selection eligibility. Tests verify root-independent
split identities and reject duplicate/missing/extra tips, negative or absent
edge lengths, and unresolved trees.

```bash
python scripts/audit_genus_codon_trees.py \
  --trees results/phylogeny/genus-codon-trees-v1 \
  --inputs results/cds/genus-codon-diagnostic-inputs-v1 \
  --output results/phylogeny/genus-codon-tree-audit-full-v1
```


The first completed-case snapshot passed all implemented audit checks for
415 cases (391 code 1, 24 code 12), 415,000 bootstrap trees and 1,686 internal
edges. Recounted UFB values agree within 0.5 percentage points. Of these edges,
806 have UFB below 95 and 468 have SH-aLRT below 80; 109 cases contain an edge
at or below 1e-5 substitutions/site. These thresholds describe uncertainty;
they are not automatic inclusion rules or calibrated biological error rates.

Warning review identifies seven cases with saturated pairwise-distance warnings,
18 with parameter-boundary warnings and 16 with unusually long NNI convergence
warnings (categories can overlap). These require review before codon-model
interpretation. Saturated nucleotide-distance warnings do not establish
synonymous saturation. Another 127 cases have repeated IQ-TREE memory-adjustment
warnings, often reporting tiny MB values; these are retained separately and
require investigation of their cause. The raw 127,301 warning lines remain in
the archived audit, outside Git. All case-level flags and summary receipts are
versioned as `metadata/genus_codon_tree_*snapshot*`.

Reproduce the categorization using `scripts/summarize_genus_tree_audit.py`
with `--audit results/phylogeny/genus-codon-tree-audit-snapshot-v1` and a new
`--output` directory. Early-finishing cases can be biased; the full audit must
be rerun after the complete queue finishes. No case, including those without
warnings, is designated selection-ready by this summary.


### Four-taxon memory-warning investigation

The previously unresolved memory-warning pattern has a source-based explanation
in IQ-TREE 3.0.1, pinned commit `d89ce077a639f5c812a10a52c022456b4e10b23a`.
In `tree/phylotree.cpp`, `LH_MIN_CONST` is 1 (line 42). In memory-saving mode,
lines 1024–1045 cap the likelihood-slot count at `leafNum - 2`, then require
at least `int(log2(leafNum) + 1)` slots. With four taxa these limits are two
and three: the warning is emitted and allocation raised to three even with
ample requested memory. `utils/tools.cpp` lines 4554–4575 confirm that `--mem 2G`
is parsed as 2 × 1,073,741,824 bytes, not two bytes or two MB.

All 127 four-taxon cases in the 415-case audit have the warning, and none of the
288 other cases does. Across 127,127 warning lines, reported adjusted allocation
is 0.011–0.894 MB. This matches the source explanation and does not indicate
exhaustion of the requested 2 GB. The inference is based on matching source and
logs; compiled-binary tracing and paired numerical reruns were not performed.
No production setting or executable was changed. Parameter-boundary, saturation
and NNI-convergence warnings still require separate review.

Source: [pinned IQ-TREE allocation implementation](https://github.com/iqtree/iqtree3/blob/d89ce077a639f5c812a10a52c022456b4e10b23a/tree/phylotree.cpp#L1024)
and [memory-option parser](https://github.com/iqtree/iqtree3/blob/d89ce077a639f5c812a10a52c022456b4e10b23a/utils/tools.cpp#L4554).
Full source files and receipt are archived in
`results/environments/iqtree-3.0.1-memory-warning-review-v1`; checksums and the
case-set comparison are versioned in `metadata/iqtree_quartet_memory_warning_review.json`.
The earlier snapshot receipt remains unchanged as a historical record.

The case-set comparison can be reproduced from the archived audit:

```python
import csv
from pathlib import Path
root = Path('results/phylogeny/genus-codon-tree-audit-snapshot-v1')
with (root / 'case_audit.tsv').open() as handle:
    quartets = {r['case_id'] for r in csv.DictReader(handle, delimiter='\t')
                if int(r['taxa']) == 4}
with (root / 'warnings.tsv').open() as handle:
    flagged = {r['case_id'] for r in csv.DictReader(handle, delimiter='\t')
               if 'Too low -mem' in r['warning']}
assert flagged == quartets and len(flagged) == 127
```


## Full MG94 diagnostic execution

The full 1,655-case information-screened queue completed in
`results/cds/genus-mg94-diagnostics-v3`. Four one-CPU workers consume completed
nucleotide-tree receipts as they become available. All 457 installed HyPhy files
are checksum-verified at startup. The producer pins the executable, upstream
`tests/hbltests/libv3/support/FitMG94.bf`, inputs, tree configuration and scripts.
The full tree/support audit remains required before interpreting codon results;
execution can proceed while that audit and the tree queue finish.

Each case uses a fixed topology with internal support labels replaced by unique
node names; root-independent splits are checked before and after export. HyPhy
2.5.101 fits a global MG94xREV model with CF3x4 frequencies, one shared omega,
and the upstream profile confidence interval. LRT is explicitly disabled.
NCBI codes 1 and 12 map to `Universal` and `Alt-Yeast-Nuclear`. The executable
passed all 128 codon-table assertions and both exact option-name assertions
against Biopython's tables; the reusable test is
`tests/hyphy_fungal_genetic_codes.bf`. This tests genetic-code interpretation,
not model adequacy.

The source exports synonymous and nonsynonymous substitution components per
site. These must not be labeled conventional dS/dN per synonymous/nonsynonymous
site without verifying the normalization. A shared-omega fit does not estimate
branch-specific omega. Complete JSON, saved likelihood-function files and logs
are retained for independent checks of tree identities, branch components,
parameter estimates and profile-interval behavior. Gene-copy and genetic-code
metadata are carried into each case configuration. Recombination, saturation,
alignment sensitivity, topology uncertainty and selection eligibility remain
unresolved by execution alone.

Planning allowances are four workers, 2 GB memory per job, 20 GB output and
2–168 hours on the existing host, with no new charges. Memory/runtime are
planning estimates rather than enforced or measured bounds. A timestamped
observation records 73 completed cases while the queue remains active. See
`metadata/genus_mg94_{resource_plan,execution_config,execution_observation}.json`.

The initial `v1` launch failed before loading the model because generic HBL
assignments require the `ENV=` argument. Its logs are preserved. The initially corrected
`v2` passed a deterministic per-case seed through `ENV=RANDOM_SEED=...;` and produced
completed fits before the later topology correction described below. No failed attempt is counted as a biological result.
The launch command is the argument interface of
`scripts/run_genus_mg94_diagnostics.py`; exact per-case commands are preserved
in each case's `config.json`. Do not launch a duplicate while active.


For a fresh run on this host, use an unused output directory. A complete tree
source is checked through its full receipt and all case hashes, without a live
producer requirement. Only an incomplete source requires --tree-producer-pid;
the original run used PID 2749871 and recorded its process start ticks.

```bash
HYPHY_SOFTWARE_ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/run_genus_mg94_diagnostics.py \
  --inputs results/cds/genus-codon-diagnostic-inputs-v1 \
  --trees results/phylogeny/genus-codon-trees-v1 \
  --information metadata/genus_codon_tree_information.tsv \
  --resources metadata/genus_mg94_preserved_topology_resource_plan.json \
  --code-check metadata/hyphy_fungal_genetic_code_check.json \
  --hyphy-install "$HYPHY_SOFTWARE_ROOT/hyphy-2.5.101-install" \
  --hyphy-source "$HYPHY_SOFTWARE_ROOT/hyphy-2.5.101" \
  --installed-manifest results/environments/hyphy-2.5.101-v1/installed_files.json \
  --output results/cds/genus-mg94-diagnostics-v3
```


### Topology-preservation correction and saved-fit readback

The v2 queue was retired cleanly after 433 completed fits. Independent readback
found that HyPhy's default `--kill-zero-lengths Yes` removed internal branches
estimated as zero during its preliminary nucleotide GTR fit: 48 of the 433
cases had fewer internal branches than their inputs. This violated the intended
fixed branch grid. All v2 outputs and the failed audit are preserved, with
case counts in `metadata/genus_mg94_v2_retirement_receipt.json`; none is mixed
into the corrected queue. Earlier execution observations are historical.

The full v3 queue explicitly passes `--kill-zero-lengths No`. Its producer now
checks both exported topology splits and the exact named branch grid before
recording case completion. The upstream option is defined in
`SelectionAnalyses/modules/shared-load-file.bf` lines 496–534 of pinned HyPhy
2.5.101. No HyPhy library or executable was changed. New resource/configuration
receipts use the `genus_mg94_preserved_topology_` prefix. The queue remains
1,655 cases; this is a full correction, not a reduced experiment.

`scripts/audit_genus_mg94_fits.py` independently reloads every saved likelihood
function in a fresh one-CPU HyPhy process and evaluates it without optimization.
It verifies case/configuration/artifact hashes, model/code options, codon
coverage, taxa, topology, named branch components, saved omega and profile
interval containment. Numerical discrepancies are retained as review flags.
The first v3 snapshot passed for 41 cases and 389 branches: maximum absolute
likelihood discrepancy 5.46e-12 and component additivity error 2e-10; no numerical
review flags were triggered. This verifies output consistency, not optimization
convergence, profile-interval calibration or biological eligibility. A full
readback remains necessary after the complete queue finishes.

```bash
python scripts/audit_genus_mg94_fits.py \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --install /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/hyphy-2.5.101-install \
  --output results/cds/genus-mg94-audit-full-v3
```


## Complete tree and MG94 execution/readback

All 1,655 nucleotide searches and all corrected v3 MG94 fits are complete.
The full nucleotide audit read 1,655,000 bootstrap trees and checked 6,721 ML
internal edges, exact taxon grids, finite branch lengths, report settings and
support frequencies. All implemented checks passed. The complete group set
contains 1,556 code-1 and 99 code-12 cases; the 57 information-screen exclusions
remain recorded separately. The full review is
`results/phylogeny/genus-codon-tree-review-full-v2`.

Tree warnings remain material review inputs: 33 cases have saturated nucleotide
pairwise-distance warnings, 53 parameter-boundary warnings, 83 NNI-convergence
warnings, and two warn that a sequence has over 50% gaps/ambiguity (Amanita and
Trichoderma, both marker 5005697at2759). Categories overlap. Memory adjustments
occur in exactly the 505 four-taxon cases, consistent with the documented
IQ-TREE slot-allocation explanation. Of 6,721 internal edges, 3,066 have UFB
below 95 and 1,723 have SH-aLRT below 80; 388 cases contain near-zero edges.
These are diagnostic summaries, not automatic biological exclusion rules.

The full MG94 readback reloaded and evaluated every saved likelihood without
optimization. All 1,655 cases and 18,407 branch records passed the implemented
input-grid, topology, parameter/interval and additive-component checks, with
zero numerical review flags. Maximum absolute likelihood discrepancy was
9.10e-12; maximum branch-component additivity error was 2e-10. This verifies
saved/reported consistency, not an optimum or calibrated uncertainty.

A separate full binding check verifies that every information-screened case
appears exactly once in both the MG94 completion grid and the independent tree
audit, and that every MG94 input retains the corresponding audited topology.
The script is `scripts/verify_genus_mg94_tree_bindings.py`; its receipt is
`metadata/genus_mg94_full_tree_binding.json`. Full audit receipts, case-level
review tables and execution summaries are versioned under the respective
`genus_codon_tree_full_*` and `genus_mg94_full_*` prefixes. Compact receipts
point to and hash their complete source receipts outside Git.

Codon-specific saturation, recombination, alignment/gene-copy sensitivity,
component normalization, optimization sensitivity and propagation of topology
uncertainty still precede selection inference. Completion of these execution
and readback stages does not establish selection eligibility.


### Pinned site-normalization definition

The upstream standalone FitMG94 analysis includes a normalization absent from
the older test-support fitter used for our completed diagnostics. It computes
synonymous and nonsynonymous expected-site opportunities S and NS, weighted by
the fitted equilibrium codon frequencies, then reports:

- dS = synonymous branch component × (S + NS) / S.
- dN = nonsynonymous branch component × (S + NS) / NS.

The [pinned upstream implementation](https://github.com/veg/hyphy-analyses/blob/42a3fd041399a17b9bc9ccb64f31f2a2fb764972/FitMG94/FitMG94.bf#L455)
uses `ComputePairwiseDifferencesAndExpectedSites(code, {})`. In the installed
HyPhy 2.5.101 helper, those defaults weight alternative nucleotide changes
equally; stop outcomes do not contribute to the S/NS numerators, while the
alternative-nucleotide denominator is retained. These are explicit upstream
counting conventions, not a claim that all possible dS definitions coincide.

Source files and hashes are archived in
`results/environments/hyphy-mg94-normalization-source-v1` and
`metadata/hyphy_mg94_normalization_source.json`. The downloaded branch tip was
verified against pinned commit 42a3fd041399a17b9bc9ccb64f31f2a2fb764972. No fitted
outputs or installed libraries were changed. The next step is to compute this
normalization from the saved fits, independently verify opportunities and
scaling under both genetic codes, then perform codon-divergence/saturation
reviews. Raw synonymous components should not yet be labeled conventional dS.


### Independent correction of stop-codon opportunity compaction

The first normalization attempt stopped before producing branch distances:
independent opportunity assertions detected a defect in the pinned HyPhy
2.5.101 `ComputePairwiseDifferencesAndExpectedSites` helper. During conversion
from 64-codon arrays to sense-codon arrays, the stop-codon branch increments the
offset but still assigns the stop's zero counts. That overwrites the preceding
sense-codon entry. In both codes 1 and 12, GTT, TAC and TCT lose their S and NS
values: 12 mismatched scalar entries across the two codes. The full comparison
is `metadata/hyphy_original_opportunity_comparison.tsv`.

The isolated correction inserts `continue` after incrementing the offset, before
assigning the S/NS values. The installed library and all saved fits remain
unchanged. Patch application with zero fuzz reproduces the isolated helper's
checksum. Every corrected opportunity entry passed independent Biopython
enumeration: 61 sense codons × two opportunity types × two codes = 244 runtime
assertions. Three regression tests cover codons preceding stops, stop-neighbor
normalization and the code-12 CTG reassignment.

The patch is `patches/hyphy-2.5.101-opportunity-compaction.patch`; the original
source hash, defect comparison and corrected runtime verification are recorded
in `metadata/hyphy_opportunity_compaction_{defect,verification}.json`. The first
attempt remains archived as `genus-mg94-normalized-branches-v1` with a failure
receipt. It is not a completed normalization or a new fitted model.

The full corrected workflow, `scripts/normalize_genus_mg94_branches.py`, reads
every saved equilibrium-codon-frequency vector, checks finite nonnegative
entries and normalization, and computes frequency-weighted S and NS using the
independent enumeration. A fresh HyPhy process for every fit loads the saved
model and calculates the same quantities with the corrected helper. Agreement
is required within 1e-12 before deriving dS and dN for that fit's branches.
The output columns explicitly name the equal-alternative counting convention.
The derived distance ratio need not equal the fitted global omega, since the
opportunity convention does not use fitted nucleotide mutation-rate weights.
It is not a separately estimated branch-specific omega.

```bash
# Apply only to a new output file; do not patch the installed library.
patch --batch --fuzz=0 \
  -o /tmp/genetic_code_corrected.bf \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/hyphy-2.5.101-install/share/hyphy/TemplateBatchFiles/libv3/tasks/genetic_code.bf \
  patches/hyphy-2.5.101-opportunity-compaction.patch
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/normalize_genus_mg94_branches.py \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --audit results/cds/genus-mg94-audit-full-v3 \
  --binding metadata/genus_mg94_full_tree_binding.json \
  --source metadata/hyphy_mg94_normalization_source.json \
  --plan metadata/genus_mg94_corrected_normalization_resource_plan.json \
  --install /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/hyphy-2.5.101-install \
  --opportunity-helper /tmp/genetic_code_corrected.bf \
  --output results/cds/genus-mg94-normalized-branches-v2
```

Use new output paths for reruns. The resource plan reserves two one-CPU workers,
4 GB memory and 1 GB output, with .1–4 hours planning time on the existing host.
Normalization does not establish synonymous-saturation limits, optimization
adequacy, time calibration, topology uncertainty or selection eligibility.


Corrected normalization is now complete for all 1,655 cases and 18,407 branches
in `results/cds/genus-mg94-normalized-branches-v2`. All 1,655 independent HyPhy
frequency-weighted opportunity checks passed; maximum disagreement with the
Python calculation was 2.23e-15. Full inverse readback verified that every
exported dS/dN component returns to its source component under inverse scaling,
and the installed original helper's checksum is unchanged. Complete provenance
and case summaries are in `metadata/genus_mg94_normalization_*` and
`metadata/genus_mg94_case_normalization.tsv`.

The median per-case maximum branch dS is 0.925 under this explicit convention,
but the largest branch estimate is 99.09. Such very large estimates require
alignment, branch-length identifiability and saturation review before biological
interpretation. These descriptive values are not calibrated eligibility
thresholds. No selection test was run, and no case is certified selection-ready
by successful normalization.


### Conditional longest-branch likelihood diagnostics

Completed seven branch-length evaluations for every case's largest normalized-dS
branch (11,585 evaluations across all 1,655 cases). Multipliers are 0.1, 0.25,
0.5, 1, 2, 4 and 10 relative to the fitted branch parameter. All other parameters
remain fixed. The script verifies actual assigned parameter values, recomputes
the baseline at multiplier 1 and restores/rechecks the original likelihood
at the end of every case. Target-node ties are resolved by node name.

The full review table records the targeted branch's terminal/internal status,
root-independent smaller-side taxon split, codon coverage, prior tree warnings
and gene-copy caveat. This connects the largest distances to their input data;
it does not label the branches as evolutionary accelerations.

No sampled point improves the saved likelihood beyond numerical precision
(maximum improvement 9.10e-12). The largest estimates are not automatically flat
along this conditional slice. For example, the Wallemia hederae terminal branch
F1540922 in marker 319730at2759 has fitted dS 99.09 and complete retained-codon
coverage. Halving, doubling and multiplying its branch parameter tenfold change
log likelihood by −53.03, −42.13 and −122.43, respectively. This does **not**
certify absence of synonymous saturation: nuisance parameters were not
reoptimized, and the global-omega model constrains synonymous and nonsynonymous
components together. Alignment/gene-copy problems can also remain despite
complete coverage and a peaked conditional curve.

![Largest-dS conditional slices](figures/genus_longest_branch_likelihood_slices.svg)

The figure shows the 12 largest estimates; every case and evaluation is retained
in `results/cds/genus-mg94-longest-branch-slices-v1`. Case reviews and compact
receipts are versioned under `metadata/genus_mg94_longest_branch_*`. These are
fixed-parameter slices, **not profile likelihoods, confidence intervals or
selection tests**. Nuisance-reoptimized profiles, alignment review and other
model sensitivities are still needed before biological eligibility decisions.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/slice_genus_mg94_longest_branches.py \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --normalized results/cds/genus-mg94-normalized-branches-v2 \
  --tree-review metadata/genus_codon_tree_full_case_review.tsv \
  --plan metadata/genus_mg94_longest_branch_slice_plan.json \
  --install /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/hyphy-2.5.101-install \
  --output results/cds/genus-mg94-longest-branch-slices-v1
python scripts/plot_genus_branch_likelihood_slices.py \
  --slices results/cds/genus-mg94-longest-branch-slices-v1 \
  --output docs/figures/genus_longest_branch_likelihood_slices.svg
```

Use new output paths for reruns. The plan reserves two one-CPU workers, 4 GB
memory, 1 GB output and .1–4 hours on the existing host.

Independent readback of all 1,655 raw HyPhy logs and 11,585 table points passed.
It verifies source hashes, deterministic target selection, multiplier grids,
branch-parameter and dS scaling, restored baselines, and summary likelihood
differences. This is a log/table audit, not a second likelihood implementation.

```bash
python scripts/readback_genus_branch_slices.py \
  --slices results/cds/genus-mg94-longest-branch-slices-v1 \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --normalized results/cds/genus-mg94-normalized-branches-v2 \
  --output metadata/genus_mg94_longest_branch_slice_readback.json
```

Choose a new readback receipt path for reruns.

## Longest-branch parameter profiles: full queue running

The next numerical diagnostic reoptimizes all remaining free MG94 parameters
while fixing the previously selected longest-dS branch's underlying `t` parameter
at 0.1, 0.25, 0.5, 1, 2, 4 and 10 times its original value. Each of the 1,655
cases also receives an unconstrained reoptimization: 13,240 fits in total.
Every point starts independently from its original saved model, uses HyPhy
`USE_LAST_RESULTS=1` and optimization precision 1e-6, and retains the same
alignment, topology, genetic code and empirical CF3x4 frequencies. Five free
exchangeabilities, global omega and the other branch parameters can change.

These are finite-grid, single-start profiles of **branch parameter t**, not dS
confidence intervals. Reoptimizing exchangeabilities changes the relationship
between t and normalized synonymous distance. Original dS remains only a review
label identifying the case and target. Optimized points may still fall below a
global profile maximum; no confidence cutoff or saturation certificate is applied.
The unconstrained fit and best evaluated point both remain explicit, so numerical
improvements over the original fit cannot be hidden by a single baseline choice.

Every optimized point is checked against its fixed-nuisance starting likelihood.
The target must stay fixed where constrained. All branch and global free values
are exported to a complete model, then a fresh HyPhy process reloads it and
recalculates the likelihood without optimization (tolerance 1e-6). The full branch
identity grid is checked against the original fit JSON before optimization.
Detailed parameters, likelihoods, elapsed times, scripts, logs and hashes are
retained outside Git; the final full-table audit remains pending.

The first launch was retired after fresh readback exposed an exporter defect:
the assignment parser omitted a last branch when a `SetParameter` command followed
its semicolon on the same line. No complete case passed that launch. The original
producer and failed artifacts remain in `genus-mg94-branch-parameter-profiles-v1`,
with a retirement receipt. The corrected producer preserves trailing commands and
checks the full branch grid. Its first completed cases pass fresh likelihood
readback; it runs in the separate immutable v2 directory. Original MG94 fits,
installed HyPhy libraries and the earlier fixed-nuisance slices are unchanged.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/profile_genus_longest_branch_parameter.py \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --slices results/cds/genus-mg94-longest-branch-slices-v1 \
  --plan metadata/genus_mg94_branch_parameter_profile_plan.json \
  --install /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/hyphy-2.5.101-install \
  --output results/cds/genus-mg94-branch-parameter-profiles-v2
```

The queue has completed; use a new output path for any later rerun. The plan
reserves four one-CPU workers, 8 GB memory, 20 GB disk and 1–24 hours on the
existing host. The original fits consumed 1.042 summed worker-hours; repeated
tighter fits and fresh readbacks justify the larger allowance. The launch receipt
is only an initial completed-case/process observation, not full-run completion.

### Source-quality caveat identified after the numerical diagnostics

The audited FCS/CDS overlap mapping identifies 29 completed genus-codon cases
containing N. cerealis (F610337) marker proteins whose coding regions overlap
publisher EXCLUDE intervals. These cases are listed in
`metadata/fcs_codon_case_exposure.tsv`. Several high-distance examples involve
this taxon, so a peaked branch likelihood does not establish a biological
acceleration signal. Source/organism identity review and explicit exclusion
sensitivities are required before biological interpretation. Existing numerical
fits and profiles remain preserved as baseline diagnostics; their likelihood
readback checks do not resolve this source-quality concern.

### FCS-omission eligibility accounted for across all 1,655 codon cases

The full-case source and taxon-grid check finds 1,626 cases unchanged by the
specified EXCLUDE/FIX/TRIM omission rule. All 29 affected Naganishia cases would
drop from four taxa to three after omitting F610337. They are therefore below the
existing four-taxon project threshold and are not refitted under that policy.
This is a design/eligibility limit, not evidence for a null effect or a general
claim that three-taxon likelihood models are impossible. Source review, better
sampling or a separately justified design would be needed for these cases.

```bash
python scripts/assess_fcs_codon_omission_eligibility.py --fits results/cds/genus-mg94-diagnostics-v3 --exposure results/qc/fcs-marker-analysis-exposure-v1 --output results/cds/fcs-omission-eligibility-v1
```

All case dispositions are versioned in `metadata/fcs_codon_omission_disposition.tsv`.
Being unchanged by this one omission rule does not establish broader selection
eligibility for the other cases.

### Completed profile audit and remaining optimization concerns

All 1,655 cases completed seven fixed-t grid points and one unconstrained
reoptimization, totaling 13,240 optimized fits and fresh saved-fit readbacks.
The full audit checked 66,200 point-artifact hashes and 226,696 fitted parameter
values. It also verified that the entire saved model text outside the profiled
parameter declarations remained unchanged, including the fixed reference
exchangeability. Maximum fresh likelihood readback discrepancy was 8.64e-11.
An additional aggregate-table check verified all 26,480 target-t and omega
values against the saved optimizer logs, with exact agreement.

Four cases retain a substantive optimization concern: a constrained grid fit
has higher likelihood than the unconstrained reoptimization by more than 1e-5.
All four are Malassezia cases (730114at2759, 4940884at2759, 541070at2759 and
4976279at2759), with excess log likelihoods approximately 1.995, 0.836, 0.538 and
0.105. These need restarts from the better solutions before using the
unconstrained optimum as a reference. They are not among the 29 FCS-exposed
cases. The numerical audit passes because it verifies what was computed;
it does not establish global optimization or selection eligibility.

Across cases, the median unconstrained likelihood improvement over the original
fit is 0.00266, with maximum 2.518. Allowing nuisance reoptimization substantially
changes some fixed-t likelihoods, reinforcing that fixed-nuisance slices are
not profile likelihoods. Neither analysis supplies a calibrated dS confidence
interval: t is the profiled parameter and exchangeabilities can change.
No saturation threshold or biological acceleration conclusion is derived here.
All 29 FCS-exposed cases are explicitly marked in the new case summary.

```bash
python scripts/audit_genus_branch_parameter_profiles.py \
  --profiles results/cds/genus-mg94-branch-parameter-profiles-v2 \
  --fits results/cds/genus-mg94-diagnostics-v3 \
  --slices results/cds/genus-mg94-longest-branch-slices-v1 \
  --plan metadata/genus_mg94_branch_parameter_profile_plan.json \
  --output results/cds/genus-mg94-branch-parameter-profile-audit-v1
```

The completed audit requires a fresh output directory on rerun. Case summaries,
receipt and target/omega table check are in
metadata/genus_mg94_branch_parameter_profile_audit_* and
metadata/genus_mg94_profile_table_field_readback.json. The full completion
receipt is metadata/genus_mg94_branch_parameter_profile_completion_receipt.json.

### Completed restarts for the four optimization discrepancies

All 32 unconstrained restarts completed, using eight previously evaluated
solutions for each flagged Malassezia case. Every saved fit underwent a fresh
likelihood replay during execution. A subsequent readback checked 194 artifact
hashes, the complete case/start grid, logged likelihoods, exported parameters,
unchanged model definitions and summary arithmetic. The maximum saved replay
discrepancy is 7.28e-12. Original profiles remain unchanged.

| Marker | Best start | Gain over previous best evaluated likelihood | Target t | Omega |
| --- | --- | ---: | ---: | ---: |
| 4940884at2759 | factor 10 | 2.255891 | 2169.302 | 0.0001211 |
| 4976279at2759 | factor 4 | 0.011193 | 230.185 | 0.0022809 |
| 541070at2759 | factor 10 | 0.179187 | 490.474 | 0.0005640 |
| 730114at2759 | factor 10 | 0.615943 | 377.675 | 0.0011325 |

Better evaluated solutions resolve the earlier ordering discrepancy but do not
establish global optima. Very large target parameters paired with small omega
require further identifiability assessment; t is not normalized dS. The finite
profile grids have not been regenerated around these solutions, and no case is
cleared for selection inference or saturation-sensitive comparisons. These four
cases are not among the separately flagged FCS-exposed cases.

Full outputs: results/cds/genus-mg94-flagged-profile-restarts-v1. Archived
receipt, points, summary and readback: metadata/genus_mg94_flagged_restart_*.
The following command checks saved artifacts without rerunning optimization;
choose a new output path for another readback:

```bash
python scripts/audit_flagged_mg94_restarts.py \
  --restarts results/cds/genus-mg94-flagged-profile-restarts-v1 \
  --plan metadata/genus_mg94_flagged_restart_plan.json \
  --output metadata/genus_mg94_flagged_restart_readback.json
```

### Full codon-case review ledger — September 23

The completed diagnostic summaries have now been joined into one explicit
review ledger covering all **1,712 prepared codon cases**: 1,655 fitted and
57 lacking sufficient sequence information for the existing tree workflow.
Case grids, taxon counts, nucleotide/codon dimensions and copy flags agree
across the source tables. The four profile-discrepancy cases map exactly to
the completed multistart refits and remain flagged for identifiability review.

There are 549 cases with at least one listed case-specific review flag and
1,163 without those listed flags. Absence of a flag is not selection eligibility.
Overlapping flag counts include 388 near-zero nucleotide-tree edges, 83 NNI
convergence warnings, 53 parameter-boundary warnings, 33 pairwise-saturation
warnings, 29 cases falling below four taxa under declared FCS omission, 14
copy-reconciliation caveats across the full input set, two other tree warnings,
and the four multistart identifiability concerns. The 57 information-limited
cases are also explicitly retained. Counts must not be summed as independent
cases, and numerical thresholds here are inherited review labels, not validated
biological exclusion rules.

Every row retains `selection_eligibility=not_established`. Saturation and
identifiability assessment, alignment/copy adequacy and the selection model/test
design remain required. This ledger supports explicit eligibility decisions
and prevents completed numerical checks from being mistaken for biological
clearance; it does not perform a selection test.

Reproduce with `python scripts/summarize_codon_analysis_readiness.py` (fresh
output directory required). Outputs are archived in
`metadata/codon_analysis_readiness.tsv` and
`metadata/codon_analysis_readiness_receipt.json`. An independent pandas join
and Boolean reconstruction reproduced the exact full case grid and flagged
case union, recorded in `metadata/codon_analysis_readiness_readback.json`.
The local summary completed in under one second without GPUs or model refits.

### Full-group realignment inputs verified — September 27

Prepared exact full source proteins and CDS for **all 1,712 existing diagnostic
groups**, covering **11,963 case–taxon sequences**. Membership, code 1/code 12
assignments, hybrid exclusions and recorded copy caveats are preserved from the
verified diagnostic inputs. The source CDS translates exactly to the original
ungapped full protein under the case's genetic code; a single terminal stop is
removed explicitly where present. Ambiguous source characters are preserved.

Every original diagnostic column is linked back to a full-protein residue index
through its recorded global MAFFT column. Independent FASTA parsing, codon-table
lookups and cumulative residue-index reconstruction passed all exports and maps,
reproducing all **5,728,643 observed original taxon–codon cells**. Original input,
source-CDS and global alignment artifact hashes were also rechecked. This
provides full-cohort inputs for alignment sensitivity beyond the earlier
20-outlier review; it does not yet perform realignment or establish selection
eligibility, correct orthology, genus monophyly or absence of saturation.

```bash
python scripts/prepare_full_codon_realignment_inputs.py --plan metadata/full_codon_realignment_input_plan_20260927.json
python scripts/readback_full_codon_realignment_inputs.py --plan metadata/full_codon_realignment_input_plan_20260927.json --output metadata/full_codon_realignment_input_readback_20260927.json
```

Use fresh output paths for reruns. Full FASTA files, original-position maps and
case summaries are under `results/cds/full-group-realignment-inputs-20260927-v1`.
The plan and completed independent proofs are archived in
`metadata/full_codon_realignment_*_20260927.json`. The input preparation allowance
was one CPU, 4 GiB RAM, 1 GiB output and 1–15 minutes; it completed locally in
seconds without GPU prediction. Whole-group alignments, codon projection and
correspondence/coverage sensitivity remain the next stages. No flagged or
information-limited case was silently excluded from the 1,712-case input set.

### All-group full-protein realignment running

Started MAFFT 7.525 `--auto --thread 1 --inputorder` for every verified group,
with four concurrent cases. This extends alignment sensitivity to the full
1,712-case set, including previously information-limited and flagged cases.
Exact input taxa/order, ungapped residues and absence of all-gap columns are
checked before case completion. Fixtures passed valid ambiguous-residue/gap
preservation and rejected duplicate/reordered taxa, substitutions, all-gap
columns and ragged outputs. Ambiguous source characters are not replaced.

The plan pins the input receipt/readback, implementation and MAFFT installation
files (59 pins including analysis inputs/scripts). Systemd enforces four CPUs,
16 GiB RAM and no swap. Pre-run planning allows 0.5–24 hours, 10 GiB output and
a 100-GiB free-disk reserve. Each subprocess has a four-hour timeout that fails
explicitly rather than silently dropping a case. This is a broad planning bound,
not a measured full-cohort ETA. Per-case commands, elapsed time, input/output
hashes and MAFFT logs are retained; completed cases are reusable only with the
same plan and verified artifacts. Unreceipted final outputs require review.

```bash
python scripts/check_full_codon_group_alignment_validation.py
python scripts/run_full_codon_group_realignments.py --plan metadata/full_codon_group_realignments_plan_20260927.json
```

The live PID/creation time and frozen plan are recorded in
`metadata/full_codon_group_realignments_{plan,launch}_20260927.json`.
Outputs are under `results/cds/full-group-protein-realignments-20260927-v1`.
Do not launch a second instance while active. These are raw alternative protein
alignments; production completion, independent readback, codon projection and
correspondence/coverage sensitivity remain pending. A change after realignment
can reflect taxon sampling, the auto-selected algorithm or site retention; it
is not evidence that either alignment is correct or a selection test is valid.

### Full realignment readback queued

Added an independent manual-FASTA audit for all 1,712 groups. It checks the
complete case universe, every source taxon/order and ungapped sequence,
rectangularity and all-gap columns, alignment dimensions, recorded code/copy
fields and commands, and all output/log hashes. The MAFFT strategy text is
recorded per case to preserve auto-selected algorithm provenance. No producer
FASTA or validation helper is imported. Valid sequence preservation and seven
malformed/corrupted cases passed fixtures.

The queued verifier checks the exact producer PID, creation time and command,
then requires its final source-bound receipt. Its allowance is one CPU, 4 GiB
RAM, no swap, 0.05 GiB output and 1–15 minutes after producer completion; systemd
enforces the resource limits. Plan and live controller records are
`metadata/full_codon_group_realignments_readback_{plan,launch}_20260927.json`.
Outputs will be under
`results/cds/full-group-protein-realignment-readback-20260927-v1`.

```bash
python scripts/check_full_codon_alignment_readback.py
python scripts/readback_full_codon_group_realignments.py --plan metadata/full_codon_group_realignments_readback_plan_20260927.json
```

The raw alignment run and full readback remain pending. Codon projection,
occupancy and correspondence comparisons require the completed audit; sequence
preservation alone does not establish homologous alignment columns or selection
eligibility.

### Full realignment verified; codon projection produced

All 1,712 realignments and their independent full readback completed successfully.
The audit covers 11,963 sequences and 6,770,822 source amino acids. MAFFT reported
L-INS-i for 1,706 cases and FFT-NS-i for six; algorithm provenance remains explicit.
Both services exited successfully, and all case/artifact hashes were rechecked.
Completed evidence is `metadata/full_codon_group_realignments_completed_20260927.json`
and its separate completed-readback file.

Projected exact source codons into every verified local alignment with recorded
code 1/code 12 assignments. Canonical translated codons count toward an integer
80% taxon-occupancy rule; ambiguous/gapped observations become whole `???`
triplets and matching `?` amino-acid states. Retained columns and full source
residue indices are exported. Fixtures cover nonstandard translation, gaps,
ambiguous codons, exact occupancy boundaries, empty masks, stops/short inputs
and a partition of correspondence losses.

The producer emitted all 1,712 cases, 41,473 taxon-pair comparisons and 6,432,298
observed taxon–codon cells. It reports all groups passing the inherited coverage
gate; this is **pending independent projection verification**, not selection
eligibility. Every original residue pair is partitioned into preserved in the
filtered local alignment, absent from the raw local alignment, or present there
but removed by local filtering. Newly retained local pairs can reflect original
site masking as well as changed correspondence; lower divergence or more columns
would not establish a better biological alignment.

```bash
python scripts/check_codon_realignment_projection.py
python scripts/project_full_codon_realignments.py --plan metadata/full_codon_realignment_projection_plan_20260927.json
```

The immutable output is
`results/cds/full-group-codon-realignment-projection-20260927-v1`. The plan and
producer summary are `metadata/full_codon_realignment_projection_*_20260927.json`.
Projection used local CPU with a pre-run one-CPU/4-GiB/2-GiB-output/1–20-minute
allowance; no model refits or GPU predictions were launched. All output hashes
passed, but the independent codon/correspondence reconstruction remains next.

### Full projected-codon and correspondence audit passed

Independent reconstruction passed every output codon, amino acid, occupancy
mask, residue-coordinate map, case summary and pair statistic: **1,712 cases,
11,963 sequences, 41,473 taxon pairs and 6,432,298 observed taxon–codon cells**.
The checker uses separate FASTA parsing, cumulative residue ranks and encoded
integer pair intersections; it does not import the producer projection helpers.
Three hundred randomized comparisons against direct sets passed before the
full audit. All cases retain the existing coverage gate, which is not selection
eligibility or a guarantee of alignment quality.

The verified tables yield these descriptive correspondence totals:

| Original residue-pair observation outcome | Count | Fraction |
|---|---:|---:|
| Preserved after local alignment and filtering | 19,791,982 | 99.0351% |
| Absent from the raw local alignment | 177,960 | 0.8905% |
| Preserved in the raw local alignment but filtered out | 14,867 | 0.0744% |
| Total original observations | 19,984,809 | 100% |

These pair observations share taxa, residues and ancestry; they are not
independent samples. Case-weighted median-of-median pair retention is 99.889%,
but the minimum case median is 36.686% (Naganishia__5013992at2759__code1).
Other low-retention cases include Tilletia__776280at2759__code1 (50.820%) and
Wallemia__5013992at2759__code1 (55.861%). At least one original correspondence
is absent from the raw local alignment in 1,155 cases; 690 have some loss due
only to local filtering. These overlapping counts are descriptive review
signals, not calibrated exclusion thresholds.

Retained columns total 917,328 locally versus 813,518 originally. Extra columns
and high pooled retention do not establish correct homology or resolve the
case-specific source-quality, copy, saturation and identifiability flags.
Low-retention cases must be joined to those existing flags rather than accepted
or rejected using a new arbitrary retention cutoff.

```bash
python scripts/check_codon_projection_readback.py
python scripts/readback_full_codon_projection.py --plan metadata/full_codon_realignment_projection_readback_plan_20260927.json
python scripts/summarize_full_codon_alignment_sensitivity.py
```

Plans, full proof, completion and descriptive summary receipts are archived in
`metadata/full_codon_realignment_projection_*_20260927.json` and
`metadata/full_codon_alignment_sensitivity_summary*_20260927.json`. The complete
case summary is outside Git at
`results/cds/full-group-alignment-sensitivity-summary-20260927-v1/cases.tsv`.
Independent pandas grouped sums checked every case count and pooled total;
grouped medians reproduced the reported median-of-case-medians. The full codon
readback used local CPU only. No selection model was fitted in this stage.
## Full realignment information screen and historical flags integrated

All 1,712 cases now have a joined ledger at
`results/cds/full-group-alignment-readiness-20260927-v1/cases.tsv`.
`scripts/integrate_codon_alignment_sensitivity.py` verifies input hashes,
reproduces every original nucleotide-information count from the original
codon FASTAs, recomputes the same screen on the local alignments, and preserves
every original readiness field alongside the full alignment-sensitivity fields.
The screen requires four distinct aligned nucleotide strings and at least one
parsimony-informative column; column counts use canonical ACGT, while distinct
strings retain missing symbols, exactly as in the original workflow.

| Original screen | Local screen | Groups |
|---|---|---:|
| Pass | Pass | 1,654 |
| Pass | Fail | 1 |
| Fail | Pass | 7 |
| Fail | Fail | 50 |

Thus 1,661 local alignments pass the information screen, versus 1,655 original
alignments. Seven changes are three-to-four distinct strings (two Microbotryum
and five Tilletia cases); Microbotryum marker 345792 changes from four to three.
All eight retain informative columns. Their full joined rows are archived in
`metadata/full_codon_alignment_information_changes_20260927.tsv`.

`scripts/readback_codon_alignment_integration.py` independently parses the
local FASTAs, recomputes all metrics with NumPy base counts, and checks every
historical and sensitivity field using exact pandas joins. The full readback
passed; receipt and proof are archived as
`metadata/full_codon_alignment_integration_receipt_20260927.json` and
`metadata/full_codon_alignment_integration_readback_20260927.json`.

All 549 historically flagged cases retain their flags. These describe the
original diagnostic fits, not newly measured problems in the local alignments.
No retention threshold was imposed and no case was cleared for selection.
Next, local tree/model sensitivity must be evaluated with explicit handling of
source-contamination and copy caveats, including the 29 FCS-omission groups
that fall below four taxa. Passing this information screen does not override
those requirements or establish homology correctness, identifiability, or
selection-test eligibility.
## Full local-alignment nucleotide trees running

The complete local-alignment sensitivity analysis now runs 1,632 supported
nucleotide-tree diagnostics. Of the full 1,712-case ledger, 51 cases fail the
existing information screen and 29 have fewer than four taxa after the declared
FCS EXCLUDE/FIX/TRIM omission policy. All cases remain in
`metadata/local_codon_tree_information_20260927.tsv`. Copy caveats and historical
fit warnings remain attached to eligible diagnostic cases; these fits do not
establish selection eligibility.

The historical FCS-disposition table covered only the 1,655 original fitted
cases. Preparation therefore checks all local FASTA taxon sets directly against
the independently audited, full marker/CDS FCS overlap table, reproduces the
original dispositions where available, and checks the 57 formerly unfitted
cases explicitly. An independent pandas/action-filter and manual FASTA join
verified all 1,712 new dispositions. Absence of a recorded overlap does not
certify a sequence as contamination-free.

`scripts/prepare_local_codon_tree_plan.py` records input, proof, executable-helper
and script hashes in `metadata/local_codon_tree_resource_plan_20260927.json`.
`scripts/run_local_codon_trees.py` is a separate variant of the existing runner;
the original remains unchanged. It uses the same deterministic per-case seeds,
GTR+F+G4 nucleotide model, retained identical sequences, 1,000 SH-aLRT replicates
and 1,000 UFB replicates with bootstrap NNI and saved trees. Four single-thread
jobs run under a four-CPU/16-GiB systemd limit with no swap. The planning allowance
is 10 GB output and 1–96 hours, not a measured completion estimate; no GPU or paid
resources are used.

Live identity is recorded in `metadata/local_codon_tree_launch_20260927.json`;
output is `results/phylogeny/local-codon-trees-20260927-v1`. Per-case receipts
verify inputs, tree tips, finite edges and bootstrap tip grids. The existing
full tree audit must subsequently verify reports, splits and support counts
before tree comparisons or downstream codon fits are interpreted. Execution,
full output audit and alignment-dependent topology/parameter comparisons are
not yet complete.
## Full local-tree audit queued

`scripts/advance_local_codon_tree_audit.py` is running behind the identified
1,632-case producer. It waits for the recorded PID/creation-time/command to
terminate, checks a complete execution receipt and the exact expected count,
then invokes the unchanged `audit_genus_codon_trees.py` without incomplete-case
mode. All scripts, configuration, input proof and disposition hashes are pinned.
An incomplete producer status and an incomplete case count were each tested
and rejected before the audit could launch.

The audit will independently parse and check all 1,632,000 saved bootstrap trees,
exact tip/split grids, reports and rounded UFB frequencies. SH-aLRT support
ranges are checked but the likelihood tests are not independently rerun.
Topology/alignment sensitivity comparisons and codon adequacy remain downstream.
The controller is limited to one CPU, 8 GiB RAM and no swap; its planning
allowance is 1 GB output and 0.1–12 hours after production, on existing resources.
Plan and live identity are recorded in
`metadata/local_codon_tree_audit_plan_20260927.json` and
`metadata/local_codon_tree_audit_launch_20260927.json`. Expected output is
`results/phylogeny/local-codon-tree-audit-20260927-v1`, with a separate controller
receipt under `results/cds/local-codon-tree-audit-handoff-20260927-v1`.
Queuing the audit does not establish that production or verification is complete.
## Alignment-dependent tree comparison prepared and queued

`scripts/compare_codon_alignment_trees.py` compares the original and local
supported nucleotide trees only after both complete tree audits pass. It
preserves the full 1,712-case ledger, with expected dispositions of 1,625
matched cases, 30 original-only cases (29 FCS exclusions and one information
failure), seven local-only cases, and 50 cases fitted under neither alignment.
Matched comparisons require exactly identical taxon sets.

For each matched case it records unrooted internal split overlap, RF distance
and RF divided by `2 * (taxa - 3)`, nucleotide tree-length changes, and all
taxon-pair path-length changes. A split-level table retains terminal and
internal edges, branch lengths, and SH-aLRT/UFB labels. Support and edge-length
deltas are computed only for shared splits; missing splits are not assigned
zero length. Inferred zero-length edges remain explicit resolutions rather
than being silently collapsed. Historical flags, copy caveats and alignment
retention accompany case summaries. Likelihoods from different alignment
columns are not compared, and these descriptive differences are not
independent samples, homology validation, or selection evidence.

`scripts/check_codon_alignment_tree_comparison.py` checks alternate orientations
of the same unrooted tree, a conflicting quartet, zero-length edges, branch
path sums against Bio.Phylo's direct distances, and invalid tip/root/length
inputs. A full original-versus-itself execution additionally checked all
1,655 audited original fits: 40,935 pair rows and 18,407 edge rows had exact zero
differences, with all splits shared. Evidence is recorded in
`metadata/codon_tree_comparison_identity_check_20260927.json`; the self-check
output is `results/cds/codon-tree-comparison-identity-check-20260927-v1`.
This validates the identity behavior, not the pending local comparison.

`scripts/advance_codon_tree_comparison.py` now waits for the identified local
audit controller and its successful full receipt before producing
`results/cds/codon-alignment-tree-comparison-20260927-v1`. The pinned plan and
live identity are `metadata/codon_tree_comparison_plan_20260927.json` and
`metadata/codon_tree_comparison_launch_20260927.json`. The controller has one
CPU/8 GiB RAM, no swap, 1 GB planned output and a 0.01–2 hour planning allowance
after auditing. No GPU or paid resources are used. The actual comparison and
its independent full output readback remain pending.
## Independent full tree-comparison readback queued

`scripts/readback_codon_tree_comparison.py` reconstructs the full comparison
with DendroPy 5.0.8, independently of the producer's Bio.Phylo parser and
split-summed distances. It checks every case disposition and retained flag,
every union split and support/length delta, and every taxon-pair distance using
DendroPy's direct phylogenetic distance matrix. Row grids and missing fields
must match exactly; floating-point comparisons use `1e-10` absolute/relative
tolerance. Input, source receipt and output artifact hashes are checked.

The checker passed the entire original-tree identity run (1,712 ledger cases,
1,655 matched cases, 18,407 split rows and 40,935 pair rows). Evidence is in
`metadata/codon_tree_comparison_identity_independent_readback_20260927.json`.
`scripts/check_codon_tree_comparison_readback.py` then confirmed that both a
wrong pair distance and a missing pair are rejected even when artifact hashes
and the declared row count are updated. No success proof was written for either
corrupted fixture.

The full realignment-comparison checker is now waiting on the recorded comparison
controller's exact PID, creation time and command. Plan and launch records are
`metadata/codon_tree_comparison_readback_plan_20260927.json` and
`metadata/codon_tree_comparison_readback_launch_20260927.json`. It uses one
CPU/8 GiB RAM, no swap, a 0.1 GB output allowance and a 0.01–2 hour planning
range after comparison. Expected proof is
`metadata/codon_tree_comparison_completed_readback_20260927.json`.
The real comparison is still pending, and the successful identity test does
not establish that the new alignments preserve topology or model adequacy.
## Full local-alignment MG94 diagnostics queued

All 1,632 information/FCS-qualified local groups are queued for MG94 refitting
after their complete independent nucleotide-tree audit. The separate runner
`scripts/run_local_mg94_diagnostics.py` preserves the original MG94xREV global
CF3x4 model, deterministic case seeds, code 1/12 handling, retained zero-length
branches and `--lrt No`. It binds the full tree audit to the exact producer
case receipts, verifies the source case metadata, and retains copy caveats.
The original runner and installed HyPhy files remain unchanged.

All 457 installed HyPhy files were checked against their manifest, and the
FitMG94 model source matches the original diagnostic fits. A fixture confirmed
that a partial tree audit is rejected before model access or fit execution.
`scripts/advance_local_mg94_diagnostics.py` waits for the identified full-tree
audit controller and its successful receipt. Plan and live identity are in
`metadata/local_mg94_handoff_plan_20260927.json` and
`metadata/local_mg94_launch_20260927.json`.

The resource plan `metadata/local_mg94_resource_plan_20260927.json` reserves
four single-thread workers, 16 GiB RAM, no swap and 20 GB output. The original
1,655 fits recorded 3,750.14 summed case-seconds (median 1.95, maximum 20.87
seconds), corresponding to 0.26 hours under ideal four-worker scheduling.
The new planning range is 0.5–8 hours after tree auditing, allowing for different
alignment lengths, numerical behavior and shared-host load; it is not a bound.
No GPU or paid resources are involved. Expected fits are under
`results/cds/local-mg94-diagnostics-20260927-v1`.

The new fits are pending, and no selection test is being run. Saved-likelihood
replay, full numerical audit, corrected post-fit site-opportunity normalization
and original-versus-local parameter comparison remain required. The previously
identified opportunity-compaction defect concerns post-fit normalization and
must use the isolated corrected helper; it does not alter the pinned MG94
fitting model. Component branch lengths must not be called conventional dS/dN
before that normalization. Source/copy, recombination, saturation and
identifiability reviews also remain open.
## Full local-MG94 likelihood replay queued

`scripts/advance_local_mg94_audit.py` waits for the identified full local-MG94
controller, verifies the exact 1,632-case completion grid and all case receipt
hashes, then runs the unchanged `audit_genus_mg94_fits.py` without incomplete
mode. The auditor reloads every saved likelihood without optimization and checks
the model/code settings, codon coverage, topology, branch identities, branch
component additivity, omega and profile interval consistency. The controller
also verifies every replay script/log hash and the exact audit-to-fit case join.

The audit preserves `complete_saved_fit_readback_with_review_flags` separately
from `passed_saved_fit_readback`; a completed controller is not an assertion
that numerical flags are absent. Neither result establishes optimality,
profile-interval calibration, biological adequacy or selection eligibility.
Corrected site-opportunity normalization remains downstream.

Plan and live identity are recorded in
`metadata/local_mg94_audit_plan_20260927.json` and
`metadata/local_mg94_audit_launch_20260927.json`. Expected output is
`results/cds/local-mg94-audit-20260927-v1`, with controller records under
`results/cds/local-mg94-audit-handoff-20260927-v1`. The queued process has one
CPU/8 GiB RAM, no swap, 2 GB planned output and a 0.1–8 hour planning range after
the fits finish. It uses existing CPU resources. Execution and numerical
validation remain pending.
## Corrected full local-MG94 normalization queued

The corrected opportunity helper formerly referenced under `/tmp` was absent.
The recorded patch was reapplied with zero fuzz to the unchanged installed
original, producing a durable isolated helper at
`results/environments/hyphy-corrected-opportunities-20260927-v1/genetic_code_corrected.bf`.
Its SHA256 exactly matches the previously verified correction
`ed0c3cb35c4e7361f751e5ccf6203bd9df098250e99354d42fce95e41f9a0aa5`.
Reproduction provenance is archived in
`metadata/local_mg94_corrected_helper_20260927.json`; no installed library changed.

`scripts/advance_local_mg94_normalization.py` now waits for the identified full
local fit audit. It requires a clean numerical audit, binds every fit to its
audited nucleotide tree with the existing full binding checker, then invokes
the unchanged normalization workflow on all 1,632 cases and 18,302 branches.
The execution plan is written only after these dependencies complete and
records their exact receipt hashes. The normalizer checks all 244 genetic-code
opportunity values and independently compares Python and corrected HyPhy
frequency-weighted opportunities for every fit. Final artifact and per-fit
check hashes are verified before the controller reports completion.

A fixture confirmed that unresolved numerical flags stop normalization before
tree binding and produce an explicit review-required record. This does not
exclude flagged cases silently or clear them for inference. The normalization
uses the same equal-alternative convention as the original fits; its derived
dN/dS ratio is not an independently fitted branch omega and may differ from
the global fitted omega.

The plan and live process are recorded in
`metadata/local_mg94_normalization_plan_20260927.json` and
`metadata/local_mg94_normalization_launch_20260927.json`. Resources are two
CPUs/4 GiB RAM, no swap, 1 GB planned output and 0.1–4 hours after auditing.
Expected results are under `results/cds/local-mg94-normalized-branches-20260927-v1`.
Normalization and original-versus-local divergence comparisons remain pending;
no GPU, paid resources, new likelihood optimization or selection test is used.
## Full local-tree execution complete; independent audit running

All 1,632 local-alignment tree fits completed, comprising 1,533 code-1 and 99
code-12 cases and 1,632,000 saved bootstrap trees. Genetic code labels remain
case provenance here; the tree fits themselves use the shared nucleotide
GTR+F+G4 model. The original 1,712-case disposition grid remains intact:
51 information-limited and 29 FCS-omission cases were not fitted.

After successful producer termination, verified the exact case grid, all 6,528
recorded tree/log/report/bootstrap artifact hashes, per-case configuration
bindings, information-table hash and pinned inputs. The completion archive is
`metadata/local_codon_tree_execution_completed_20260927.json` (producer receipt
SHA256 `f1a5ba9c56653fdd08061252fc83243302125c13ce9c5ddacc0d1230a8d42949`).
Summed recorded case execution time was 8,069.55 seconds; this is not wall time
for the whole pipeline.

The existing independent audit controller automatically advanced and was
verified live by PID, creation time and command. Its full report, split and
support reconstruction is not yet complete. The queued tree comparison and
codon refits retain their audit dependency; this production checkpoint does
not establish topology stability, model adequacy or selection eligibility.
## Full normalized-divergence sensitivity comparison queued

`scripts/compare_codon_alignment_divergence.py` now waits for completed local
normalization, then compares the audited original baseline and local fits under
the same equal-alternative site-opportunity convention. It verifies full clean
fit audits, exact fit/normalization case bindings, artifact hashes, code labels
and tree-node grids. Matched cases must have identical taxon sets and genetic
codes. All 1,712 ledger cases remain visible, including unmatched/unfitted cases.

The output records case-level global omega and total normalized dS/dN tree
distances, shared-branch distance differences using canonical unrooted splits,
and all taxon-pair path-distance differences. Internal node labels are never
used to match branches across fits. Absent splits remain missing rather than
receiving zero length, and their differences are not computed. Path distances
can be compared across topology changes on the exact same taxon pairs.

A complete original-versus-itself check passed for all 1,655 fitted cases,
18,407 branch rows and 40,935 pair rows, with every difference exactly zero.
The plan and evidence are `metadata/codon_divergence_identity_plan_20260927.json`
and `metadata/codon_divergence_identity_readback_20260927.json`. This checks
identity behavior, not the pending realignment results; independent full
comparison reconstruction remains required.

The actual comparison plan and live identity are
`metadata/codon_alignment_divergence_comparison_plan_20260927.json` and
`metadata/codon_alignment_divergence_comparison_launch_20260927.json`.
Output is planned under `results/cds/codon-alignment-divergence-comparison-20260927-v1`,
using one CPU/8 GiB RAM, no swap, 1 GB output and 0.01–2 hours after normalization.
No GPU or paid resources are used. Original baseline fits are not substituted
with later multistart optima; historical warnings remain attached explicitly.
These are descriptive, correlated sensitivity measurements, not cross-alignment
likelihood tests, branch omega estimates or evidence of positive selection.
## Local-tree audit and alignment sensitivity completed

The full independent local-tree audit passed all 1,632 fits, 1,632,000 saved
bootstrap trees and 6,703 internal edges. Model reports, exact tip/split grids
and rounded UFB frequencies were checked; SH-aLRT ranges were checked without
rerunning the tests. The completed receipt is archived in
`metadata/local_codon_tree_audit_completed_20260927.json`. Numerous warning log
lines remain in the audit output and are not equivalent to that many unique
problematic cases or a declaration of model adequacy.

The original-versus-local comparison also completed and passed full independent
DendroPy reconstruction: 1,712 ledger cases, 18,705 union split rows and 40,755
taxon-pair rows. There are 1,625 matched cases, 30 original-only, seven local-only
and 50 fitted under neither alignment. Completion is recorded in
`metadata/codon_tree_comparison_completed_20260927.json` and
`metadata/codon_tree_comparison_completed_readback_20260927.json`.

**315/1,625 matched groups (19.38%) change inferred unrooted topology.** Among
changed groups, normalized RF has median 0.25 (range 0.142857–1). There are
448 original-only and 448 local-only internal splits; 6,243 internal splits and
all 11,566 terminal splits are shared. Inferred zero-length resolutions remain
in these counts.

Original-only internal splits have median UFB 49% (quartiles 41–58%, range
2–99%); local-only splits have median 53% (quartiles 44–65%, range 10–100%).
Most changed splits therefore have modest support, but support is not uniformly
low. These marginal distributions do not by themselves identify strongly
supported incompatible split pairs; that requires explicit incompatibility
checks. No statistical significance is inferred from these correlated branches.

Of the changed groups, 149 have historical review flags and 166 do not. Median
case-level residue-pair retention is 0.997573 among changed cases and 1.0 among
unchanged cases. High alignment retention therefore does not guarantee identical
inferred topology, and absence of a historical flag does not establish
robustness or eligibility. This comparison also changes retained alignment
columns, so it is not a controlled equal-column likelihood comparison.

`scripts/summarize_codon_tree_sensitivity.py` produces the full changed-case
table, support-frequency table and historical-flag cross-tabulation under
`results/cds/codon-tree-sensitivity-summary-20260927-v1`. Receipt and independent
count/quantile checks are archived as
`metadata/codon_tree_sensitivity_summary_20260927.json` and
`metadata/codon_tree_sensitivity_summary_readback_20260927.json`.
The local MG94 controller advanced automatically and had 70 completed fits at
the checkpoint. Divergence sensitivity and selection eligibility remain open.

## September 27 paired split incompatibility

The completed original/local comparison was checked for direct incompatibility
rather than comparing marginal support maxima. All 804 old-only/new-only split
pairs across the 315 topology-changing cases were examined. Exactly 589 pairs
are incompatible, spanning all 315 cases; none of these pairs contains a
zero-length branch. Four nonempty intersections define incompatibility for
unrooted splits. A separate checker independently enumerated resolved quartets
for every candidate pair and verified the entire output grid, quartet witnesses,
branch supports and lengths, historical flags, alignment retention, and support
curve.

The maximum, across incompatible pairs, of the smaller of the two ultrafast
bootstrap supports is **85%**. Thus the earlier 99%/100% marginal maxima do not
represent a conflict supported at those levels on both trees. This does not
establish that either alignment or topology is correct, nor remove optimization,
recombination, model-adequacy or selection-eligibility concerns. Multiple pairs
within a case are correlated. The entire empirical support curve is retained;
no threshold was used to select the analyzed cases.

Reproduce from the checksum-bound comparison with:

```bash
python scripts/summarize_codon_split_conflicts.py --output results/cds/codon-split-conflicts-20260927-v1
python scripts/readback_codon_split_conflicts.py
```

The producer requires a new output directory. Detailed split pairs are in that
results directory. Provenance, complete support curve and independent readback
are archived in `metadata/codon_split_conflicts_*_20260927.*`.

## Independent normalized-divergence comparison readback

`scripts/readback_codon_alignment_divergence.py` reconstructs branch identities
with DendroPy, replaces tree edge lengths with each audited normalized dS/dN
metric, and calculates direct phylogenetic distance matrices. It checks every
case, split and taxon-pair row, global omega, tree totals, missing cells,
historical flags, translation code, and exact fit/audit/normalization grids.
Source and output checksums are checked before and after reconstruction.
Numerical comparisons use absolute and relative tolerances of 1e-10.

The complete baseline identity comparison passed independently: 1,712 ledger
rows, 1,655 fitted groups, 18,407 splits and 40,935 taxon pairs. The proof is
`metadata/codon_divergence_identity_independent_readback_20260927.json`.
The corruption checks in `scripts/check_codon_divergence_readback.py` alter
pair distance, pair completeness, and global-omega difference while updating
artifact/plan hashes, ensuring semantic checking rather than hash-only rejection.

Once the actual comparison producer finishes, run:

```bash
python scripts/readback_codon_alignment_divergence.py --plan metadata/codon_alignment_divergence_comparison_plan_20260927.json --proof metadata/codon_alignment_divergence_completed_readback_20260927.json
```

The actual local/original comparison remains pending. This checker validates
comparison arithmetic against audited fitted values, not optimization quality,
model adequacy, saturation clearance or selection eligibility.

The full independent comparison readback is now queued automatically through
`scripts/advance_codon_divergence_readback.py`. The controller follows the
recorded producer PID, creation time and command, waits for termination, then
requires a complete checksum-bound comparison receipt before running the
checker. Its pinned plan is
`metadata/codon_divergence_readback_handoff_plan_20260927.json`; its live launch
identity is recorded separately. Limits are one CPU, 8 GiB RAM, no swap, and
10 MB expected output, with a planning allowance of 0.01–2 hours after the
producer. The terminal-producer path was exercised on the full baseline identity
dataset and passed. Actual results will be written under
`results/cds/codon-divergence-readback-handoff-20260927-v1`; the controller fails
if any required source or the expected full disposition grid is missing.

## September 27 full local MG94 execution and numerical audit completed

All 1,632 local-alignment MG94 fits completed (1,533 genetic-code-1 and 99
code-12 groups). The full independent saved-likelihood audit passed across
18,302 branches with zero numerical consistency flags. Maximum absolute
reported-versus-replayed log-likelihood difference is 7.275957614183426e-12;
maximum branch-component additivity error is 2e-10. Summed per-case fitting
time was 3,783.35 seconds, distinct from wall time under four workers.

After both execution and audit services exited successfully, rechecked the
complete case/configuration grid, all 6,528 fit-artifact hashes, all 3,264
likelihood-replay script/log hashes, both audit tables, genetic-code counts,
and explicit `not_established` selection eligibility for every case. Completion
proof is `metadata/local_mg94_completed_audit_20260927.json`.

Fit receipt: `07b5fd7405ab6d7e9a7e602f118273cd3f9ef0097e59053eb5a1cbff570ea439`.
Audit receipt: `cc272632d553d9ac7be625edcecfc42530b3470939e9d00e9a2ccaf34f43a2aa`.
These results establish saved-fit numerical consistency only. They do not prove
optimality, calibrated profile intervals, model adequacy, absence of saturation,
or suitability for selection testing. Corrected opportunity normalization and
the full original/local comparison remain downstream.

## Local diagnostic review ledger

`results/cds/local-codon-diagnostic-review-20260927-v1/cases.tsv` now joins the
complete local tree/fit diagnostics to all 1,712 original ledger cases. All
historical flags and original columns are retained; only the historical-scope
annotation is updated. A separate dataframe reconstruction checks every column
and every review flag. Scripts are `integrate_local_codon_diagnostics.py` and
`readback_local_codon_diagnostics.py`; receipts are archived under
`metadata/local_codon_diagnostic_review_*_20260927.json`.

Among local cases, 365 have near-zero nucleotide-tree edges under the existing
<=1e-5 review label, 593 have tree warning lines, 96 have fit warning lines,
14 retain copy caveats, and 80 have no local fit. These categories overlap;
934 cases have at least one listed local review trigger. There are no local
>=10 nucleotide-tree-edge labels or saved-fit numerical-consistency flags.
Warning lines can repeat and do not count independent problems. Fit warning
lines do not contradict the successful numerical replay audit. Their scientific
implications still require review. The broad tree-warning category here is not
directly comparable to the earlier categorized historical flag count.

All 549 historically flagged cases retain those flags, including earlier
optimization concerns. Absence of current listed flags does not clear those
concerns or establish biological adequacy. Selection eligibility remains
`not_established` for every case.

## September 27 normalized divergence sensitivity completed

Corrected opportunity normalization completed for all 1,632 local fits and
18,302 branches. The original/local comparison and its independent DendroPy
readback both exited successfully. The full comparison includes 1,712 ledger
cases: 1,625 matched, 30 original-only, seven local-only and 50 neither-fitted.
Every one of 18,705 split rows and 40,755 taxon-pair rows was independently
reconstructed. Proof is archived as
`metadata/codon_alignment_divergence_completed_readback_20260927.json`.

The case-level sensitivity summary is generated by
`scripts/summarize_codon_divergence_sensitivity.py` under
`results/cds/codon-divergence-sensitivity-summary-20260927-v1`. An independent
csv/math/statistics reconstruction checked every source field, ratio, zero
case, within-case pair median, topology flag and all 15 quantile rows. The
summary receipt, quantiles and readback are archived in metadata.

Across matched groups, median log2(local/original) is 0.003241 for total
normalized dS and 0.220850 for total normalized dN: back-transformed ratios
1.00225 and 1.16542, respectively (approximately +0.2% and +16.5%). Global
omega has median log2 ratio 0.216039 across 1,624 positive/positive groups
(back-transformed 1.16154); one zero/zero case remains explicitly undefined
for log ratios. No pseudocount is introduced. Pair-distance sensitivity is
summarized within each group before calculating across-group quantiles.

These shifts combine changes in alignment coverage, residue correspondence,
topology and fitted parameters. They do not identify which alignment is correct
or demonstrate positive selection. Groups share species and protein families;
these are descriptive distributions, not independent replicates or significance
tests. Original fits remain the archived baseline, including their unresolved
optimization concerns. All historical warnings and eligibility limitations
remain in force.

## Reproducible alignment-sensitivity figure

![Codon alignment sensitivity](figures/codon_alignment_sensitivity_20260927.png)

[PDF](figures/codon_alignment_sensitivity_20260927.pdf),
[SVG](figures/codon_alignment_sensitivity_20260927.svg), and
[data](figures/codon_alignment_sensitivity_20260927.tsv) are available.
Panel A shows full empirical distributions of log2 local/original tree-total
dS, tree-total dN and fitted global omega. Panel B includes every matched case,
showing dN sensitivity versus median retention of original residue-pair
correspondences, colored by topology change. One zero/zero omega ratio remains
undefined; no pseudocount, outlier clipping or causal trend fit is used.

Reproduce with `scripts/plot_codon_alignment_sensitivity.py` (new artifacts only).
`scripts/readback_codon_alignment_figure.py` recomputes every plotted value from
the raw audited comparison and checks all case identities, flags and topology
labels. The PNG was visually inspected for full ranges, axis and legend labels,
counts and caveats. The figure does not establish either alignment as correct,
selection, or independent replication across gene groups.

## September 27 local fit warning content resolved

All 96 warning blocks from the 1,632 local MG94 fit logs are HyPhy advisories
about identical aligned sequences. The full audit in
`scripts/audit_local_codon_duplicate_warnings.py` checks each fit log and input
alignment against its recorded hashes, matches every warning to this advisory,
and independently reconstructs the reported redundant-sequence count from
aligned FASTA strings. It found 96 identical-sequence groups, containing 196
taxon memberships and 100 redundant aligned sequence copies. The complete
1,632-case table is under
`results/cds/local-codon-duplicate-warning-audit-20260927-v1`; the summary and
identical-group inventory are archived in metadata.

These advisories are not optimizer-failure messages. This content classification
does not establish optimization success or remove historical optimization
concerns. Identical strings are relative to the retained alignment, which can
include missing characters; they do not prove full protein or genome identity.
All species remain in the diagnostic analysis and no near-zero-edge, copy,
independence or selection-eligibility concern is cleared. The earlier generic
fit-warning flags remain traceable in the review ledger.

## Local longest-branch identifiability diagnostics started

The realigned fits now undergo the same longest-branch diagnostic sequence as
the original fits. The first stage evaluates the saved branch-time parameter
at multipliers 0.1, 0.25, 0.5, 1, 2, 4 and 10 for every one of the 1,632 fitted
cases, selecting the branch with greatest normalized dS in each case. All
nuisance parameters remain fixed at this stage. This gives 11,424 conditional
likelihood evaluations plus 1,632 restored-baseline checks. The branch identity,
codon coverage, genetic code, copy caveat and local tree warning count remain
recorded. No cases are selected by an arbitrary dS cutoff.

A separate script, `scripts/slice_local_mg94_longest_branches.py`, reuses the
original slice procedure while retaining the local audit's raw warning count
instead of historical warning categories. Original scripts/results remain
unchanged. Source fit, normalization and tree-review hashes are bound in
`metadata/local_mg94_longest_branch_slice_plan_20260927.json`; verified process
identity is archived in the matching launch record. Output is
`results/cds/local-mg94-longest-branch-slices-20260927-v1`. Limits are four CPUs,
8 GiB RAM and no swap, with 2 GiB output and a 0.1–4-hour planning allowance.
Existing CPU resources incur no new charges.

These fixed-nuisance slices are a prerequisite for optimized nuisance-parameter
profiles. They do not establish an optimum, provide calibrated dS confidence
intervals, clear saturation, or establish selection eligibility. Full slice
verification and optimized profile follow-up remain required.

## September 27 full local branch-profile optimization

Verified all 1,632 completed local longest-branch slices, including all
11,424 raw likelihood points, longest-normalized-dS target choices, parameter
scaling, source hashes and restored baselines. Maximum conditional likelihood
improvement is 7.28e-12 (numerical rounding). Proof:
`metadata/local_mg94_longest_branch_slice_readback_20260927.json`.

Launched all 1,632 cases for eight independent optimizations each: one
unconstrained refit and seven fixed branch-t values, with remaining branch
parameters, exchangeabilities and omega reoptimized. Each of the 13,056 fits
gets a fresh saved-model likelihood readback. Four workers, four CPUs, 8 GiB
and no swap; planning allowance 1–24 hours and 20 GiB output. Verified the
corrected service and its worker children live; the first launch failed before
output creation on an incorrect installation path, and its journal is retained.
Plan and identity are in `metadata/local_mg94_branch_parameter_profile_{plan,launch}_20260927.json`.

The finite, single-start t grid is not a normalized-dS confidence interval or
a selection result. Full profile audit remains pending; the historical auditor
has a hardcoded old-cohort FCS exposure source and must be adapted with verified
local case provenance before use. All eight aims remain open; GPU predictions
remain paused.

## September 27 local fitted-case FCS exposure verified

`scripts/audit_local_codon_fcs_exposure.py` binds all 1,632 local fits to their
hashed configuration, alignment and saved model. All 11,599 case–taxon rows
and 18,846,954 embedded nucleotide/missing characters match the source FASTA.
The full 59,840-row audited marker map joins every fitted entry; zero retained
entries overlap reported EXCLUDE/FIX/TRIM actions. This is evidence about the
existing FCS calls, not absence of contamination or selection eligibility.
Full rows and proof are in `results/qc/local-codon-fcs-exposure-20260927-v1`;
the receipt is archived as `metadata/local_codon_fcs_exposure_readback_20260927.json`.

A separate `scripts/audit_local_branch_parameter_profiles.py` now accepts an
explicit `--exposure` directory and checks the full case universe and positive
subset. It preserves the historical profile auditor and running producer.
The new audit has passed CLI parsing but awaits the complete profile batch;
its full numeric/artifact validation is not yet run or queued. The producer
was revalidated by PID, creation time and command with four worker children.

The local profile validation is now automatically queued via
`scripts/advance_local_branch_profile_audit.py`, with pinned scripts, source
receipts, exposure tables and producer identity. It requires successful terminal
service state and the entire 13,056-fit batch before auditing. Outputs will be
`results/cds/local-mg94-branch-profile-audit-20260927-v1` and the corresponding
`local-mg94-profile-audit-handoff-20260927-v1` controller record. At launch the
service was verified waiting; completed profile validation remains pending.
