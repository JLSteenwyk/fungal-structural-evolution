# Ancestral reconstruction case inputs

Status: ancestral sequences and structures have not been inferred. This stage
identifies candidate nodes for the 13 whole/domain contrast cases; it does not
replace the requirement to reconstruct ancestral ensembles for selected families.

All 26 resolved gene trees (13 families under each of the profile and MAFFT
species guides) were extracted from the fully checked reconciliation outputs.
For each case we retained the focal duplicate MRCA and its first three ancestral
levels. All 104 neighborhoods, 3,722 descendant records and 52 guide comparisons
passed independent parent-path reconstruction and exact source-tree checks.
All 52 comparisons have identical descendant sets across guides. The two
reconciliations share underlying gene-tree information, so agreement is a
sensitivity result rather than independent support.

The largest neighborhood examined per case is shown below. These are candidate
sampling neighborhoods, not selected final alignment sets or reconstructed nodes.

| Family | Case taxon | Full family proteins | Level-three proteins | Level-three taxa |
| --- | --- | ---: | ---: | ---: |
| OG0000054 | Heliocybe sulcata | 5775 | 11 | 3 |
| OG0000095 | Rhodonia placenta | 4610 | 83 | 24 |
| OG0000107 | Cryoendolithus antarcticus | 3458 | 9 | 8 |
| OG0000152 | Scedosporium apiospermum | 4020 | 7 | 4 |
| OG0000230 | Furculomyces boomerangus | 2215 | 24 | 3 |
| OG0000294 | Jaapia argillacea | 1992 | 9 | 8 |
| OG0000336 | Synchytrium microbalum | 1735 | 10 | 8 |
| OG0000972 | Neolecta irregularis | 1548 | 622 | 389 |
| OG0001082 | Leucosporidium creatinivorum | 940 | 111 | 96 |
| OG0001200 | Smittium simulii | 762 | 29 | 9 |
| OG0001203 | Leucosporidium creatinivorum | 1297 | 21 | 13 |
| OG0002650 | Piloderma croceum | 842 | 77 | 24 |
| OG0002812 | Phycomyces blakesleeanus | 1695 | 12 | 11 |

Node labels remain guide-specific; comparisons use descendant identities. A
matched ancestral level alone would not establish node homology if descendants
differed. All focal MRCAs contain exactly the two nominated duplicate proteins.
No bootstrap support, rooting confidence or branch timing follows from this
membership check.

Next required work is to recover the exact extant sequences, assess full-family
and local alignment quality and domain architecture, and choose supported
reconstruction nodes with appropriate outside sequences. Ancestral inference
must assess model, alignment and topology sensitivity and propagate sitewise
uncertainty into sequence ensembles. Indel uncertainty needs separate treatment.
Selecting only reference-agreeing cases would bias the study; experimental
reference disagreements and missing coverage remain explicit selection factors.
GPU structure prediction remains paused.

Reproduce with `python scripts/prepare_case_ancestral_neighborhoods.py`, then
`python scripts/readback_case_ancestral_neighborhoods.py`. Both preserve existing
outputs. Local source trees, descendant lists and tables are in
`results/ancestral/case-neighborhoods-20260927-v1/`; source and artifact checksums
are in its receipt. Verified closure:
`metadata/case_ancestral_neighborhoods_completed_20260927.json`.

## Exact extant sequences recovered

Recovered all 1,025 candidate-clade proteins across 417 taxa and 13 families
(482,007 residues). Each exact string matched both the checksum-bound staged
reconciliation FASTA and its independently identified QC proteome record. A
separate export readback checked all expected family/gene identities and every
sequence hash and length. Duplicate sequences and all copies are retained.

| Family | Proteins | Minimum length | Maximum length | Unique sequences | Noncanonical proteins |
| --- | ---: | ---: | ---: | ---: | ---: |
| OG0000054 | 11 | 218 | 248 | 11 | 0 |
| OG0000095 | 83 | 213 | 749 | 82 | 0 |
| OG0000107 | 9 | 942 | 975 | 9 | 0 |
| OG0000152 | 7 | 159 | 496 | 7 | 0 |
| OG0000230 | 24 | 220 | 762 | 19 | 0 |
| OG0000294 | 9 | 460 | 1091 | 9 | 0 |
| OG0000336 | 10 | 245 | 433 | 10 | 0 |
| OG0000972 | 622 | 50 | 1009 | 616 | 1 |
| OG0001082 | 111 | 205 | 921 | 111 | 0 |
| OG0001200 | 29 | 148 | 595 | 28 | 0 |
| OG0001203 | 21 | 110 | 308 | 21 | 0 |
| OG0002650 | 77 | 55 | 483 | 77 | 0 |
| OG0002812 | 12 | 770 | 817 | 12 | 0 |

Noncanonical record: `F157183_CAD6936802.1` in OG0000972, symbols/counts `{"X": 3}`. The original sequence is preserved; its alignment and ancestral-state treatment must be explicit.

Inputs are in `results/ancestral/case-sequences-20260927-v1/`. Reproduce with
`python scripts/recover_ancestral_case_sequences.py`; it refuses to overwrite
existing outputs. Closure: `metadata/case_ancestral_sequences_completed_20260927.json`.
This recovery used one CPU, streamed existing local FASTAs, and required no GPU
or network access. Alignment, outside-clade context and ancestral ensembles
remain pending.

## Two-method alignment run started

All 13 candidate sets are queued through MAFFT 7.525 L-INS-i (`--amino
--localpair --maxiterate 1000 --thread 4 --threadit 0`) and FAMSA 2.2.3
(`-t 4 -keep-duplicates`). Runs are sequential with four threads each, an 8 GiB
service memory limit and no swap. The 26-alignment planning envelope is
0.25–24 hours and 1 GiB output; no GPU or paid resources are used.

Every output must contain exactly the input identifiers and reproduce every
original sequence after gap removal, including the three X residues. No copies
are removed, including identical sequences. Column occupancies at 50%, 70% and
90% are descriptive exports; no trimming is applied and occupancy is not a
measure of alignment homology correctness. The first family passed preservation
checks under both methods at launch verification; remaining runs are pending.

`scripts/run_ancestral_case_alignments.py` uses the frozen
`metadata/ancestral_case_alignment_plan_20260927.json`; exact process identity
is in `metadata/ancestral_case_alignment_launch_20260927.json`. Completed runs
can be reused only with matching settings and checksums. An unfinished run
directory is preserved and requires explicit recovery rather than silent
overwriting. Full independent output audit, residue-correspondence comparison,
domain/fragment assessment and reconstruction sampling decisions remain pending.

## Full alignment readback and correspondence comparison queued

The downstream checker waits for successful terminal completion of all 26
alignments, validating the producer identity while it runs. It will independently
reconstruct every original sequence and residue position, every coverage column
and all threshold counts. All IDs, source receipts and artifact checksums must
match the frozen input and execution records.

For each family, exact full-column correspondence between MAFFT and FAMSA uses
the complete vector of residue positions (or gaps) across all proteins. This
criterion is deliberately stringent: one changed protein can break a full-column
match. All columns, occupancies and matched alternative column numbers are
retained. Focal duplicate residue-pair correspondences are compared separately
so that full-clade and focal-pair sensitivity can be distinguished. Neither
agreement fraction is a probability that a column is homologous. No trimming
is applied based on these diagnostics.

The checker uses one CPU and 4 GiB RAM, no swap, with a 0.5 GiB output allowance
and a planned 0.1–2 hours after alignment completion. Script:
`scripts/audit_ancestral_case_alignments.py`. Plan and exact live identity:
`metadata/ancestral_case_alignment_audit_plan_20260927.json` and
`metadata/ancestral_case_alignment_audit_launch_20260927.json`.
Full validation and alignment-sensitivity results remain pending.
