# Robust domain candidate annotations — September 27, 2026

All 947 cross-guide stable event/domain candidates at the 30-residue/70%-coverage
screen and descriptive 0.1 Å margin now have Pfam 38.2 names, descriptions,
types and clan identifiers. This retains all 739 concordant, 202 discordant and
six sequence-unresolved cases, including outgroups. No functional class or
quality threshold was used to choose which candidates to annotate.

For each candidate, exact gene/node identities resolve all tied references,
annotation policies and domain boundaries. Duplicate interval triads shared
between guides are counted once. The 1,777 distinct triads contribute 56,864
fitting alternatives (32 per triad). Every included alternative passes the
original n30/c70 screen and unique-fit requirement. The annotation table reports
minima and maxima for common-core size, original-domain coverage, duplicate-pair
RMSD and identity, reference-distance contrast, each protein's mean core pLDDT,
and the fraction of shared residue triples with all three pLDDTs at least 70.
Triad identifiers are retained for traceable follow-up.

## Descriptive findings

There are 194 fungal discordant event/domain combinations plus eight outgroup
combinations. The most represented Pfam annotations among the fungal cases are:

| Pfam | Annotation | Combinations | Families | Taxa |
| --- | --- | ---: | ---: | ---: |
| PF00067.28 | Cytochrome P450 | 7 | 2 | 6 |
| PF00106.32 | Short chain dehydrogenase | 4 | 3 | 4 |
| PF00505.26 | HMG box | 4 | 2 | 4 |
| PF00069.32 | Protein kinase domain | 3 | 2 | 3 |
| PF00082.28 | Subtilase family | 3 | 1 | 3 |

These are candidate counts, not enrichment or independent evolutionary origins.
Pfam descriptions are family annotations, not experimentally established
functions of these particular proteins. All 583 role/class/Pfam summary groups
remain in the output, including concordant and outgroup groups.

Across the 194 fungal discordant combinations, the minimum absolute structural
contrast over all alternatives ranges from 0.103 to 1.305 Å (median 0.170 Å).
Minimum shared-core length ranges from 30 to 489 residues (median 138).
The minimum fraction of shared triples with all three pLDDTs at least 70 ranges
from 0.729 to 1 (median 0.970). Every candidate's minimum mean core pLDDT for
all three proteins exceeds 70, but this does not calibrate structural error.
Minimum duplicate-pair core sequence identity spans 0.292 to 1 (median 0.752).
Identical aligned cores still require review of flanking context, full-chain
sequence and prediction differences; they do not establish sequence-independent
structural evolution.

## Reproduction and verification

Run `scripts/annotate_robust_domain_candidates.py --plan
metadata/duplication_domain_candidate_annotation_plan_20260927.json`, then
`scripts/check_robust_domain_candidate_annotations.py --plan
metadata/duplication_domain_candidate_annotation_plan_20260927.json --output
<new-readback-path.json>`. Preserve the immutable output directory
`results/structural_comparisons/duplication-domain-candidate-annotations-20260927-v1`.
It contains `candidate_annotations.tsv`, `pfam_summary.tsv` and a hash-bound
receipt; large tables remain outside Git.

Independent dataframe joins and aggregation checked every copied source field,
Pfam annotation, triad membership, metric extremum and summary denominator.
The annotation file also matches the recorded Pfam 38.2 release checksum.
The initial checker exposed a trailing space in one Pfam description; its
independent parser now trims surrounding whitespace consistently, and the full
readback passed. No scientific values were changed. See
[completion evidence](../metadata/duplication_domain_candidate_annotation_completed_20260927.json).

The margin is a descriptive screen, not a biological significance or error
threshold. Contrasts describe similarity to a reference, not ancestral change
or additive branch distances. Case review, prediction-source controls,
nonduplication backgrounds and phylogenetically controlled tests remain necessary
before mechanistic or duplication-associated conclusions.

## Subsequent correspondence sensitivity

The [sequence-locked shared-reference control](duplication-sequence-locked-reference-controls-20260927.md)
now tests all 48 candidates with identical complete domain intervals. Four no
longer maintain the original direction beyond the descriptive margin across
all alternatives: three span both signs and one enters the margin band. The
947 candidates above remain the historical original-control set, not a claim
that all pass every subsequent control. Original fields are preserved alongside
new-control status in the subsequent candidate summary.
