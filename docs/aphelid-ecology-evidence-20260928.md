# Assembly-linked aphelid ecological evidence

The versioned ecological table now contains 47 taxa, preserving all 45 prior
rows and adding two submitter-reported, sample-linked annotations for a sparsely
sampled lineage. Both selected genomes have an explicit assembly-report link to
the BioSample, BioProject and named strain listed below.

| Selected taxon | Assembly | Strain | BioSample | Reported role and host |
|---|---|---|---|---|
| Amoeboaphelidium occidentale | GCA_023515795.1 | FD01 | SAMN26806364 | Algal parasitoid; Scenedesmus dimorphus |
| Amoeboaphelidium protococcarum | GCA_023515745.1 | FD95 | SAMN26794594 | Algal parasitoid; Scenedesmus dimorphus |

The [FD01 BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN26806364) and
[FD95 BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN26794594) describe
isolation from algal cultures in New Mexico and cultivation with the named host.
The corresponding [FD01 BioProject](https://www.ncbi.nlm.nih.gov/bioproject/817714)
and [FD95 BioProject](https://www.ncbi.nlm.nih.gov/bioproject/817550) corroborate
the reported role. These pages contain submitter annotations, not an independent
reanalysis of the experiments. The selected-isolate experimentally-verified
flag therefore remains false, while the assembly-to-sample linkage is recorded
explicitly in its own table. Neither taxon is counted as an independent origin
of parasitism, and neither is silently coded as a negative ECM observation.

Raw assembly reports and HTML snapshots are stored outside Git in
`data/traits/aphelid-ecology-20260928`. Exact URLs and SHA256 values are in
`config/aphelid_ecology_evidence_20260928.json`. Browser retrieval initially
failed for the NCBI pages, but direct HTTPS downloads succeeded and the full
saved descriptions were inspected. No ecological claim depends on a search
snippet. The source-linked genome article remains a literature lead; its full
text was unavailable during this review and is not claimed as inspected.

Run `python scripts/extend_aphelid_ecology_evidence.py` to reproduce the new
version using the pinned snapshots (existing outputs are protected). It checks
selected identities, assembly accessions, exact sample/project/strain links,
source descriptions, and preservation of every existing row. Outputs:

- `metadata/species_ecology_evidence_aphelids_20260928.tsv`
- `metadata/aphelid_ecology_sample_links_20260928.tsv`
- `metadata/aphelid_ecology_evidence_receipt_20260928.json`

Existing 32- and 45-taxon phylogenetic and structural-coverage analyses remain
bound to their original inputs. Coverage assessment and explicit trait coding
for this expanded set remain to be performed; this addition does not complete
the ecological association aim.

## Full expanded structural coverage completed

Both source-specific assessments now cover all 47 taxa and all 1,081 unordered
pairs. Independent FASTA parsing and dense mask multiplication verified every
taxon, pair, pair-marker and provisional-group count. Every original 45-taxon
row was recovered exactly in the expanded tables.

A. occidentale has 109 eligible ESMFold markers with 31,841 qualified paired
AA/3Di observations; A. protococcarum has none in either qualified source
collection. Neither aphelid has eligible AlphaFold markers in these frozen
inputs. These are alignment-eligibility findings, not claims that public
predictions or recovered protein sequences do not exist.

Across all 47 taxa, AlphaFold covers 34 taxa and 545 pairs with a marker sharing
at least 50 qualified columns; ESMFold covers 33 taxa and 528 pairs. In the
bookkeeping union, 46 taxa have some qualified coverage and 878 pairs meet the
within-source criterion; 195 pairs qualify in both and 203 in neither. No
cross-source observations were pooled. The extra 32 eligible ESMFold pairs do
not represent independent parasitism transitions or an ecological effect.

Reproduce preparation with `scripts/prepare_aphelid_ecology_coverage.py`, run the
existing assessment/readback scripts using
`metadata/qualified_{afdb,esmfold}_aphelids47_ecology_plan_20260928.json`, then run
`scripts/summarize_aphelid_ecology_coverage.py`. Full tables remain under
`results/ecology/qualified-{afdb,esmfold}-aphelids47-overlap-20260928-v1`.
The completion receipt is `metadata/aphelids47_ecology_coverage_completed_20260928.json`.
Each source retained the one-CPU, 4-GiB, 1–10-minute planning allowance and
finished in seconds without GPU inference.

## Why A. protococcarum lacks qualified single-copy markers

A source-bound trace identifies an upstream copy-selection limitation. Its
125-marker raw BUSCO table contains one Complete marker, 119 Duplicated markers
and five Missing markers. All 119 duplicated markers map to multiple annotated
gene identifiers in the completed annotation audit. The frozen single-copy
marker mapping therefore contains only marker 5001734at2759, protein
KAI3647604.1. This is not evidence that the other119 marker families lack
protein sequences, and duplicated hits alone do not distinguish biological
duplication, hybrid ancestry or assembly redundancy.

The selected protein already has an ESMFold model. In the frozen paired
alignment,53 of177 columns are qualified; the existing eligibility rule requires
max(50,ceil(0.3*177))=54. The other columns comprise13 noncanonical/missing
sequence positions, one invalid native feature and110 low-feature-confidence
positions. There is no missing structural mapping or high-PAE exclusion in
this particular row. No confidence or eligibility thresholds were changed.

Thus another round of prediction on the current single-copy list would not
address the main119-marker exclusion. The next step is to inventory all of
those copies, retain annotated gene identities and assess them through
copy-aware family phylogenies and structural comparisons. Selecting an
arbitrary best-scoring paralog would erase the uncertainty that caused the
exclusion. The existing whole-proteome family analyses are the appropriate
place to connect this information.

Reproduce the trace with `scripts/diagnose_aphelid_marker_coverage.py`; raw-table,
annotation and qualified-mask source hashes and exact counts are recorded in
`metadata/aphelid_marker_coverage_cause_20260928.json`.

## All-copy inventory completed

The source-bound inventory retains all241 BUSCO hit proteins from120 detected
marker families, including all119 duplicated markers. There are235 distinct
amino-acid sequences; the241 protein identities remain separate even where
sequences are identical. Every hit has an assignment in each of the two
existing family partitions:482 assignment rows, zero missing assignments,
and117 distinct family labels per guide. BUSCO marker labels and these family
labels are different classifications and their counts should not be equated.

`inventory_aphelid_marker_copies.py` verifies the raw BUSCO hit universe,
annotation decisions, source sequences, full family-database hash and completed
family-database readback. It exports all sequences plus per-protein gene IDs
and both family assignments under
`results/ecology/aphelid-marker-copy-inventory-20260928-v1`. FASTA readback
preserved every source sequence and protein identity. No arbitrary best-copy
selection or new prediction was performed. Structure availability, family-tree
placement and explanation of the multiple copies remain downstream work.

During implementation, a receipt-field mismatch was corrected before output
creation. An initial unindexed per-protein database lookup was stopped before
outputs and replaced by queries using the existing guide/gene index. No source
database was modified. The final inventory completed successfully; its hashes
and counts are in `metadata/aphelid_marker_copy_inventory_20260928.json`.

### Copy placements in retained gene trees

`assess_aphelid_copy_tree_placements.py` verifies inventory artifact hashes,
the retained-tree catalog checksum, and every used tree checksum. All241
proteins occur in117 validated trees, covering121 marker/family groups per
guide. The two partitions agree for all28,920 unordered pairs of target
proteins; that comparison concerns this target subset, not all family members.
The two guide rows use the same retained trees and are not independent evidence.

Per guide,117 multi-copy marker/family groups form an exclusive unrooted split,
and two do not:5003022at2759 inOG0001748 and5013067at2759 inOG0000653.
Marker530740at2759 has three copies split betweenOG0000510 (two copies,
exclusive split) andOG0001807 (one copy). The other singleton group is the
previously identified single-copy marker5001734at2759. All copies are retained.
These are descriptive placements without bootstrap qualification, duplication
age, hybrid-origin inference or arbitrary ortholog selection. Within-family
patristic distances use source-tree branch units, not structural displacement.

Outputs are in `results/ecology/aphelid-copy-tree-placements-20260928-v2`;
`metadata/aphelid_copy_tree_placements_20260928.json` records counts and hashes.
The v1 exploratory output is preserved; v2 adds explicit catalog-receipt
checksum enforcement and reproduces identical output-table hashes. Every
serialized table was read back and compared with its computed records.
