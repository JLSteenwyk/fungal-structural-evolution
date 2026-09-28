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
