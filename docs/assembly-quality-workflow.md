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

## Coding-region overlaps and affected analyses

The full mapping and independent audit are complete for all 86 assemblies with
validated flagged regions. The audit enumerated every interval against the
4,479,554 annotated CDS feature rows in those assemblies, finding 28,150
region/CDS-feature intersections involving 11,092 distinct taxon–protein IDs.
Of those proteins, 11,062 are retained gene representatives. Counts include
review-only flags and are not confirmed contaminant-gene counts. Intersections
use CDS segments, not complete gene spans across introns. Per-protein overlap
lengths union overlapping segments separately by FCS action and sequence.
All flagged sequence IDs were represented in the GFF files in this dataset.

The original 59,840 marker–taxon observations include 54 with coding overlap:

| Taxon | Marker observations | FCS action |
|---|---:|---|
| Naganishia cerealis (F610337) | 44 | EXCLUDE |
| Acaulospora colombiana (F27376) | 9 | EXCLUDE |
| Amphimedon queenslandica (O400682) | 1 | REVIEW |

The N. cerealis report assigns the affected sequences to the FCS budding-yeast
division. This is a consequential source-quality hypothesis, not independent
confirmation that those proteins are contaminants. Some large codon-distance
estimates involve this taxon; biological acceleration interpretation must await
source review and exclusion sensitivities. Taxonomic misidentification and
report/database errors remain possible alternatives to mixed-source assembly.

Exact protein-ID and sequence-hash joins locate affected observations in frozen
structural snapshots, and hashed FASTAs establish their actual use:

- Earlier ESMFold paired alignments: 16 N. cerealis marker observations.
- Combined ESMFold paired alignments: the same 16 plus three A. colombiana
  observations, 19 in total.
- Expanded AlphaFold paired alignments: none of these 54 observations present.
- Completed genus MG94 inputs: 29 cases include affected N. cerealis markers.

The first ESMFold sequence–structure coupling results therefore also require
sensitivity analysis for these observations. The absence of this particular
exposure in the AlphaFold inputs does not establish absence of other quality
problems. Baseline results and running inputs remain unchanged. These flags will
support explicit comparison with sensitivity datasets and should not be silently
converted into whole-species or whole-protein-family exclusions.

The independent checker verifies all positive and negative region/CDS
intersections using array comparisons, annotation protein IDs, interval union
lengths, representative decisions and all marker-overlap joins. Source GFF hashes
and embedded assembly-version headers are checked. Boundary checks cover inclusive
endpoints, adjacent non-overlap, intron-only overlap and duplicate segment unions.
This stage does not yet map the exact affected coding bases onto structural
residues, confirm foreign origin, or refit affected phylogenies.

```bash
python scripts/map_fcs_regions_to_cds.py --fcs results/qc/selected-assembly-fcs-audit-v1 --annotations metadata/annotation_download_receipts.json --representatives metadata/gene_representatives_receipt.json --markers results/phylogeny/markers-full-v1/protein_mapping.tsv --marker-index results/cds/marker-source-index-v1 --output results/qc/fcs-cds-overlap-v1
python scripts/audit_fcs_cds_overlap.py --mapping results/qc/fcs-cds-overlap-v1 --fcs results/qc/selected-assembly-fcs-audit-v1 --annotations metadata/annotation_download_receipts.json --representatives metadata/gene_representatives_receipt.json --markers results/phylogeny/markers-full-v1/protein_mapping.tsv --output results/qc/fcs-cds-overlap-audit-v1
python scripts/trace_fcs_marker_analysis_exposure.py --mapping results/qc/fcs-cds-overlap-v1 --audit results/qc/fcs-cds-overlap-audit-v1 --fits results/cds/genus-mg94-diagnostics-v3 --output results/qc/fcs-marker-analysis-exposure-v1
```

Full intersection and review tables remain in these result directories; compact
receipts, the 86-taxon mapping summary, 54 affected marker observations and
analysis-exposure tables are versioned under `metadata/fcs_*`.

## Existing alternative marker hits reviewed

To assess whether affected marker coverage can be recovered, reviewed every
existing raw HMMER hit and final BUSCO call for all 54 flagged marker–taxon
observations. The review verifies original protein lengths and sequence hashes,
representative decisions, source BUSCO-table hashes, HMM query identities and
coordinates, and complete HMMER report trailers. It does not rerun searches or
replace any selected sequence.

There are 1,004 raw marker–protein hits in these reports. No alternative final
BUSCO call lacks a recorded CDS/FCS overlap. Raw outputs contain 567 alternative
hits without recorded overlap; 41 of those meet or exceed their marker's absolute
dataset score cutoff, covering 35 affected marker–taxon observations. The 41
candidates comprise 32 N. cerealis and nine A. colombiana proteins. Their HMM
profile union coverage ranges from 0.2411 to 0.9987 (median 0.8056); coverage and
copy identity therefore need individual review.

The installed BUSCO 6.1.0 source also contains a relative-score filter: within a
hit category, matches below 85% of the best score can be removed. Source hashes
and the precise function are recorded in
`metadata/busco_relative_hit_filter_source_review.json`. This demonstrates why
absolute cutoff passage is not equivalent to final acceptance; it does not
reconstruct the individual rejection reasons for these candidates.

All 41 exact candidate protein sequences were materialized for review and fully
read back against their original hashes. Existing strict CDS audits report 37
exact translations and four non-triplet CDS records. These checks do not establish
orthology, lack of contamination or a suitable replacement. Full gene-copy and
phylogenetic placement review, alignment coverage, and CDS qualification remain
necessary. No candidate has been substituted into a baseline or sensitivity run.

```bash
python scripts/review_fcs_marker_alternative_hits.py --output results/qc/fcs-marker-alternative-hits-v1
python scripts/prepare_fcs_alternative_review_sequences.py --review results/qc/fcs-marker-alternative-hits-v1 --output results/qc/fcs-alternative-review-sequences-v1
```

Raw-hit tables and candidate FASTA remain outside Git; compact receipts, all 54
marker summaries and the 41 candidate review rows are versioned under
`metadata/fcs_*alternative*` and `metadata/fcs_alternative_candidate_review.tsv`.

### Exploratory marker-specific phylogenetic discordance

A completed tree for marker 129234at2759 (472 taxa, 166 retained profile-alignment
columns) places the FCS-flagged N. cerealis protein KAJ9108408.1 on a five-taxon
unrooted split with Candida albicans, C. dubliniensis, C. maltosa and C. tropicalis.
The separating branch has reported SH-aLRT support 87.2 from 1,000 replicates;
this is not bootstrap support or a probability of contamination. Tree/input
hashes, tip grids, dimensions and finite branches were verified, and the split
was identified explicitly from both sides of each unrooted edge.

This contrasts with the concatenated guides' nearest Naganishia neighbors and
supports investigating individual marker provenance. It does not by itself
distinguish contamination, misannotation, paralogy, other gene histories or model
error, and it does not justify reassignment of the whole assembly. This case was
selected for exploratory source review after the FCS/high-distance findings.

```bash
python scripts/review_flagged_naganishia_marker_tree.py --output results/qc/naganishia-flagged-marker-tree-review-v1
```

The versioned receipt is `metadata/naganishia_flagged_marker_tree_review.json`.
Other marker placements and proposed alternatives remain to be evaluated.

### Alternative candidate profile alignment and independent readback

All 41 candidate proteins were aligned to the same pinned HMMER marker models
used in the baseline profile matrix and projected onto its unchanged retained
sites. They span 33 markers and 35 affected marker–taxon observations. Their
retained-site coverage ranges from 19.97% to 100%, with median 74.63%; this differs
from raw HMM-domain coverage because the phylogenetic matrix retains a subset of
profile match states. No candidate is entirely missing.

An independent parser read the Stockholm alignment blocks, reconstructed original
protein coordinates, and checked all 14,798 candidate site rows (9,855 canonical
residues). It also verified all 5,757,596 baseline characters in the augmented
alignments. Each augmented alignment preserves all 526 baseline entries and adds
candidate gene copies with separate identifiers. Across the 33 alignments, 1,465
baseline rows have no observed amino acids; the audit records these explicitly
for removal before tree inference. These are per-marker rows, not 1,465 taxa.

Candidates remain hypotheses requiring gene-copy and placement review. Short
coverage can weaken placement, and the four non-triplet CDS records remain
unqualified for codon analyses. No replacement or new phylogeny was accepted at
this stage. Candidate sequences are not additional species.

```bash
python scripts/align_fcs_alternative_candidates.py --candidates results/qc/fcs-alternative-review-sequences-v1 --matrix results/phylogeny/profile-matrix-50-v1 --profiles results/phylogeny/profile-alignments-full-v1 --output results/qc/fcs-alternative-profile-alignments-v1
python scripts/audit_fcs_candidate_alignments.py --source results/qc/fcs-alternative-profile-alignments-v1 --candidates results/qc/fcs-alternative-review-sequences-v1 --matrix results/phylogeny/profile-matrix-50-v1 --output results/qc/fcs-alternative-profile-alignment-audit-v1
```

Compact receipts and all 41 coverage rows are versioned in
`metadata/fcs_alternative_profile_alignment_*`. Full alignments, residue mappings,
and the 17,399-row sequence-disposition table remain outside Git.

### Exploratory alternative-copy trees launched

Applied the established individual-marker coverage rule to both baseline proteins
and alternative copies: at least `max(50, ceil(0.30 * retained columns))` canonical
residues. Of 41 candidates, 38 qualify across all 33 markers. Three candidates
remain excluded from this tree run, rather than being discarded from the review:
F27376/CAG8452709.1 (5003856at2759), and F610337/KAJ9110660.1 and
KAJ9113858.1 (5011777at2759). Their observed fractions are 24.20%, 19.97% and
26.07%, respectively. Original selected proteins remain in the review alignments
when they pass the same coverage rule, preserving the ability to compare copies.

Input verification covered all 17,399 source sequence rows and 5,284,204 retained
characters. It excludes 1,519 sequence rows across these marker alignments,
including all entirely missing rows and the three below-threshold candidates.
The full disposition table is stored with the inputs outside Git.

The 33 exploratory trees use IQ-TREE 3.0.1 model selection among LG/WAG/JTT with
empirical frequencies and Gamma rates, 1,000 SH-aLRT replicates, and preservation
of identical tips. Two concurrent jobs use two threads and a 4 GB memory limit
each. Existing host memory and disk headroom were checked before launch. The
resource plan uses 48 completed baseline marker-tree runtimes: median-based
projection about 20.8 wall hours, planning range 8–72 hours, with 20 GiB disk
headroom. This is a runtime estimate, not a guarantee or new paid resource.

```bash
python scripts/prepare_fcs_candidate_gene_trees.py --source results/qc/fcs-alternative-profile-alignments-v1 --audit results/qc/fcs-alternative-profile-alignment-audit-v1 --output results/phylogeny/fcs-candidate-copy-tree-inputs-v1
python scripts/run_fcs_candidate_copy_trees.py --inputs results/phylogeny/fcs-candidate-copy-tree-inputs-v1 --output results/phylogeny/fcs-candidate-copy-trees-v1 --resources metadata/fcs_candidate_copy_tree_resources.json
```

Launch identity, input receipts, marker dispositions and resource estimates are
versioned under `metadata/fcs_candidate_copy_tree_*`. Tree completion, model and
support audits, uncertainty in candidate placement, broader paralog sampling,
and orthology decisions remain pending. A supported placement alone will not
establish an acceptable replacement or resolve assembly provenance.

### Accessibility covariates for FCS omission sensitivity

Prepared exact baseline solvent-accessibility rows for all retained observations
in the changed sensitivity markers, preserving unchanged alignment columns,
sequence characters and masks. The earlier ESMFold cohort has 223,145 observed
residues across 16 markers, 1,738 marker–taxon cells and 1,725 models; the combined
cohort has 427,510 residues across 18 markers, 3,056 cells and 3,014 models.
These cohorts overlap and must not be treated as independent datasets.

Full readback checked every retained row against the original raw accessibility,
confidence, atom count and context fields, residue mappings and paired AA/3Di
alignments. Both existing reference scales were then applied and independently
checked for every normalized row, preserving all original fields and rebuilding
scale/threshold totals. No terminal residue required withholding normalization.
At the diagnostic relative-accessibility threshold 0.25, the scales classify
9,248 earlier-cohort rows (4.14%) and 17,408 combined-cohort rows (4.07%)
differently. This measures normalization-convention sensitivity; these counts
are not independent biological changes or tests of evolutionary association.

The full command manifest is `metadata/fcs_accessibility_execution.json`, with
source receipts, full readbacks and compact summaries under
`metadata/{esmfold,esmfold_combined}_fcs_accessibility_*`. The new subset script
requires exact identity of retained columns and characters and stops if that
condition changes; it does not recompute solvent area. Existing coordinate
quality and missing-partner limitations still apply.

These tables cover the 16/18 changed markers. The 56/71 unchanged markers reuse
baseline accessibility. Later site-level exposure summaries must be recalculated
across the retained sensitivity taxa for affected markers before coupling models
are refit. Sensitivity rate estimates and revised coupling tests remain pending.
