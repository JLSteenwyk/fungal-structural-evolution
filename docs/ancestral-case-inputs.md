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
