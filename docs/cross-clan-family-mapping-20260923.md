# Cross-clan candidate family mapping

The full inventory of 68 boundary/cluster/Pfam-pair candidates now has source
protein, taxon, and family assignments under both the profile and MAFFT guides.
The 4,600 candidate-member rows expand to 4,746 protein-linked rows, representing
2,759 distinct source proteins across 224 taxa. All members are retained, before
confidence or direct-alignment screening; boundary views overlap.

Results are outside Git in
`results/structural_clusters/cross-clan-family-mapping-20260923-v1/`:

- `candidate_proteins.tsv`: complete member-to-protein identities and both family assignments.
- `candidate_family_summary.tsv`: 136 rows covering all 68 candidates, separately
  for all members and pair/cluster-exclusive models. Counts include each side
  and the intersection for models, proteins, taxa, and both family assignments.
- `receipt.json`: source hashes, output hashes, scope, and counts; archived as
  `metadata/cross_clan_family_mapping_completed_receipt_20260923.json`.

Among exclusive-model summaries, two candidate entries share a family between
the two Pfam sides under each guide. This is descriptive assignment overlap,
not proof of homology, orthology, duplication, or a structural transition.
Subsequent case selection must join the direct-comparison confidence results,
inspect domain architectures, and evaluate the corresponding gene trees.

Validation independently rebuilt the complete export using a SQL join and
compared all 4,746 serialized identities. Pandas/set aggregation checked every
count in all 136 summary rows. An initial audit query was stopped because its
query plan scanned the candidate table repeatedly; explicit join order completed
the same full validation. No analysis producer was interrupted.

Reproduce with one CPU (observed runtime seconds, excluding the inefficient
initial audit), no GPU, and local input databases:

```bash
python scripts/map_cross_clan_families.py \
  --candidates results/structural_clusters/cross-clan-candidate-inventory-20260923-v1 \
  --composition results/structural_clusters/domain-cluster-composition-20260923-v1 \
  --composition-audit metadata/domain_cluster_composition_completed_readback_20260923.json \
  --output results/structural_clusters/cross-clan-family-mapping-20260923-v1
```

The output directory must not already exist.

## Direct-comparison join within shared families

`scripts/screen_cross_clan_shared_families.py` joins every exclusive-model
cross-clan pair within an assigned family to the completed bidirectional
structural comparisons. The two qualifying candidate entries are the alignment
and envelope views of PF00085.27 versus PF26973.1 in the same structural cluster.
Both guides assign the shared family as OG0000030.

The complete join contains 1,312 protein-pair/guide/boundary rows representing
631 distinct interval pairs. These repeated views and source-protein links are
not independent observations. For each guide separately, rows passing minimum
TM-score 0.5 and minimum coverage 0.8 are:

| Boundary | No confidence fraction cutoff | 0.5 | 0.8 | 0.9 |
|---|---:|---:|---:|---:|
| Alignment | 172 | 172 | 172 | 169 |
| Envelope | 375 | 371 | 367 | 342 |

Confidence fractions require pLDDT at least 70 in both whole domains and matched
residues in both alignment directions. No PAE qualification is included. These
screens prioritize gene-tree and domain-architecture inspection; they are not
tests of evolutionary transitions or homology.

All serialized pair identities were independently reconstructed with pandas
merges, all direct scores matched their source table, and every screen flag was
recomputed. Results reside in
`results/structural_comparisons/cross-clan-shared-family-screen-20260923-v1/`;
the source/output hash receipt is archived in
`metadata/cross_clan_shared_family_screen_completed_receipt_20260923.json`.
The join uses one CPU, takes seconds, and requires no new structural inference.

```bash
python scripts/screen_cross_clan_shared_families.py \
  --mapping results/structural_clusters/cross-clan-family-mapping-20260923-v1 \
  --comparisons results/structural_comparisons/cross-clan-confidence-summary-20260923-v1 \
  --output results/structural_comparisons/cross-clan-shared-family-screen-20260923-v1
```

## Annotation and taxonomic review (September 24)

The local Pfam 38.2 annotation calls PF00085 "Thioredoxin" and PF26973
"DJC16, Thioredoxin-like domain 2", despite their different clan assignments.
Thus this candidate does not establish a previously unknown fold relationship.
It remains useful for checking annotation consistency and family history.

The PF26973 side comprises *Rozella allomycis* protein EPZ34216.1
(`24_3474`, model `AF-A0A075AV72-F1-v6`) in both boundary views, with
*Pneumocystis jirovecii* XP_018228094.1 (`365_324`, model
`AF-A0A0W4ZDZ6-F1-v6`) additionally represented in the envelope view.
Both are ingroup taxa in the sampling table. Numerous pairwise comparisons
therefore do not provide numerous independent changes on this side of the
comparison. No replicated transition or branch acceleration is established.

The exact three boundary/taxon/protein/model records, Pfam annotation rows,
and hashes of the complete source tables are recorded in
`metadata/cross_clan_shared_family_annotation_review_20260924.json`.
Full domain architecture and gene-tree inspection remain pending.

## Focal gene-tree context (September 24)

Membership/context inspection now places all 283 candidate-left proteins and
both candidate-right proteins in the 5,118-tip OG0000030 resolved tree under
each guide. The right-side proteins have a common ancestor containing 500 tips
in each native rooted tree; this count does not establish a domain transition.
Their immediate sister labels are identical under both guides:

- *Rozella* EPZ34216.1: `F2606893_XP_031856394.1`.
- *Pneumocystis* XP_018228094.1: `F263815_XP_007871995.1` and
  `F4754_XP_018226490.1`.

None of those immediate sister proteins is among the candidate-left proteins
in this cross-clan screen. The next annotation review therefore needs these
sister proteins as well as the focal proteins, without assuming that an absent
screen membership means an absent domain. Both guides giving the same local
context is a sensitivity result, not independent replication.

`scripts/inspect_shared_family_gene_trees.py` performs this inspection with
Biopython, checks both complete source-tree file hashes against the prior
membership audits, and records all inputs in
`metadata/cross_clan_shared_family_tree_context_20260924.json`. An independent
DendroPy parser confirmed both full tree sizes, focal branch lengths, immediate
sister memberships, parent labels, and the 500-tip common-ancestor scope.
An attempted ETE check was unavailable because that package was not installed;
no environment change was needed. These are checks of existing tree context,
not new inference, support assessment, rooting validation, reconciliation-event
validation, or evidence of structural acceleration.

## Full focal and sister annotation context (September 24)

All five proteins have full-proteome candidate annotations. The four overlap
policies agree on retained hits for each protein:

| Protein | Ordered retained Pfam annotations |
|---|---|
| F281847_EPZ34216.1 | PF00085.27, PF26973.1 |
| F2606893_XP_031856394.1 | PF00085.27 |
| F42068_XP_018228094.1 | PF00085.27, PF26973.1, PF24541.2 |
| F263815_XP_007871995.1 | PF00085.27, PF24541.2 |
| F4754_XP_018226490.1 | PF00085.27 |

Both focal proteins therefore contain both Pfam models in the original
cross-clan comparison. The three immediate sister proteins have PF00085
annotations despite not belonging to the candidate-left set. The original
screen restricts membership by structural cluster and model exclusivity;
its membership cannot be substituted for whole-protein domain content.
The missing retained PF26973 calls in the sisters do not yet establish domain
absence or gains/losses. Raw-hit review, sequence completeness, and comparisons
of corresponding regions remain necessary before interpreting a transition.

The complete annotations, coordinates, policy alternatives, sequence identifiers,
and input hashes are archived in
`metadata/cross_clan_focal_architectures_20260924.json`. All five exports were
checked using an independent two-step protein/sequence lookup, and the database
hash matched the earlier full-proteome snapshot. Reproduce with one CPU and no
GPU (database hashing dominates this short lookup):

```bash
python scripts/inspect_focal_domain_architectures.py \
  --context metadata/cross_clan_shared_family_tree_context_20260924.json \
  --architectures results/domains/full-candidate-architectures-v1/candidate_architectures.sqlite \
  --output metadata/cross_clan_focal_architectures_20260924.json
```

The output file must not already exist. This export records candidate annotation
context; it does not validate biological domain architectures or infer events.

## Precompetition hit review

The five proteins have 12 stored hits before overlap competition. Nine are
retained; the three excluded hits are PF13728.13 (TraF), overlapping the stronger
PF00085 hits in XP_007871995.1, XP_018228094.1, and XP_018226490.1. Their annotation
description does not establish plasmid-transfer function in these proteins.
No stored PF26973 hit was removed from any of the three sister proteins.
Thus overlap filtering does not explain their missing retained PF26973 calls;
search sensitivity and the corresponding sequence regions still need review.
These stored hits already passed the original search reporting thresholds.

Protein lengths in these hits are 510/512 residues for the Rozella focal/sister
pair and 382/387/384 for the Pneumocystis focal/two sisters. Similar lengths alone
do not validate completeness or correspondence of individual domains.

`metadata/cross_clan_focal_raw_hits_20260924.json` preserves all 12 complete hit
records, retained/excluded identifiers, and source hashes. Every retained
annotation field was checked against its raw hit, and protein/sequence links and
raw-hit counts matched the preceding export. Reproduce with:

```bash
python scripts/inspect_focal_raw_domain_hits.py \
  --architectures metadata/cross_clan_focal_architectures_20260924.json \
  --database results/domains/full-domain-database-v1/domains.sqlite \
  --output metadata/cross_clan_focal_raw_hits_20260924.json
```

## Source sequences for region comparisons

All five full-length source sequences were extracted into
`results/domains/cross-clan-focal-sequences-20260924-v1/proteins.faa`.
Each sequence exactly matches its domain-annotation SHA-256 identifier;
independent FASTA readback also confirmed every raw-hit protein length.
The archived receipt is
`metadata/cross_clan_focal_sequences_receipt_20260924.json`.
This prepares the corresponding-region comparison; it is not an alignment or
a completeness assessment. Reproduce with:

```bash
python scripts/prepare_focal_protein_sequences.py \
  --architectures metadata/cross_clan_focal_architectures_20260924.json \
  --proteomes data/qc_proteomes \
  --output results/domains/cross-clan-focal-sequences-20260924-v1
```

## Corresponding-region alignment sensitivity

MAFFT 7.525 local-pair and global-pair iterative alignments of all five proteins
were run separately with one thread and 1,000 maximum iterations. This small
CPU-only check took less than a second including hashing; no prediction was run.
Input proteins were recovered exactly by removing alignment gaps.

For the focal PF26973 alignment boundaries, the Rozella comparison maps 85 of
91 focal residues to sister residues (30 identical under local-pair, 29 under
global-pair). Each Pneumocystis comparison maps 72 of 73 residues, with 32 and
28 identical residues respectively. Envelope boundaries map 88/99 and 95/96
residues respectively. All Pneumocystis residue correspondences are identical
between settings. The Rozella mappings share 84/85 alignment-boundary pairs
and 87/88 envelope-boundary pairs.

Thus these alignments contain corresponding sequence across most of the focal
domain span in the sisters without retained PF26973 annotations. This motivates
structural and profile-sensitivity checks rather than a domain-loss call.
Repeated thioredoxin-like domains can still be misaligned, and agreement of two
settings does not independently validate homology or ancestral domain content.

All 12 explicit residue maps and identity counts were independently recalculated
from aligned prefixes. Source/software hashes, exact commands, output hashes,
and all mappings are archived in
`metadata/cross_clan_focal_alignment_receipt_20260924.json`; alignments and logs
reside in `results/domains/cross-clan-focal-alignments-20260924-v1/`.

```bash
python scripts/align_focal_domain_regions.py \
  --sequences results/domains/cross-clan-focal-sequences-20260924-v1/proteins.faa \
  --architectures metadata/cross_clan_focal_architectures_20260924.json \
  --context metadata/cross_clan_shared_family_tree_context_20260924.json \
  --mafft-root /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions \
  --output results/domains/cross-clan-focal-alignments-20260924-v1
```
