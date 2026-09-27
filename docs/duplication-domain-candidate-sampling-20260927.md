# Domain candidate sampling — September 27, 2026

All 77,760 cross-guide comparison rows have been joined to the frozen,
reconciliation-eligible taxonomy for 526 taxa. The summaries distinguish
event/domain combinations, distinct family/gene pairs, families, Pfams and taxa.
The 2,916 lineage rows include explicit zeros for all 27 role/lineage groups,
six screens, three margins and six candidate classes. The 36,941 family rows
cover observed groups only. Broad lineage groups come from the manifest;
they are not newly inferred clades or independent evolutionary origins.

At the 30-residue/70%-coverage screen and descriptive 0.1 Å margin:

| Relationship | Event/domain combinations | Gene pairs | Families | Pfams | Taxa | Fungal/outgroup combinations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Concordant in both guides | 739 | 698 | 433 | 412 | 97 | 715/24 |
| Discordant in both guides | 202 | 199 | 159 | 152 | 67 | 194/8 |
| Stable structure, unresolved sequence | 6 | 6 | 6 | 6 | 4 | 6/0 |

Distinct counts can overlap between classes and must not be summed as independent
observations. Among the 194 fungal discordant combinations, 84 are assigned to
Basidiomycota and 51 to Ascomycota. This is descriptive concentration, not
enrichment: model coverage, family composition and taxon sampling differ greatly.
Reconciled taxon counts are availability denominators, not ascertainment weights.
The contrasts measure relative similarity to a reference, not ancestral change,
absolute rates or statistically established decoupling. The 0.1 Å margin is not
a calibrated biological or prediction-error threshold.

## Reproduction and verification

Run `scripts/summarize_domain_candidate_sampling.py --plan
metadata/duplication_domain_candidate_sampling_plan_20260927.json`, then
`scripts/check_domain_candidate_sampling.py --plan
metadata/duplication_domain_candidate_sampling_plan_20260927.json --output
<new-readback-path.json>`. Output creation is exclusive; preserve existing results.

Outputs are under
`results/structural_comparisons/duplication-domain-candidate-sampling-20260927-v1`:
`annotated_comparisons.tsv`, `lineage_summary.tsv`, `family_summary.tsv` and
`receipt.json`. Large tables remain outside Git; the plan pins input locations
and hashes, and the producer receipt hashes all output tables.

The independent checker verified every source field, taxonomy assignment,
classification and exact aggregation set, including all zero-count lineage rows.
See [completion evidence](../metadata/duplication_domain_candidate_sampling_completed_20260927.json).
Candidate annotation, nonduplication backgrounds and phylogenetically controlled
tests remain outstanding.
