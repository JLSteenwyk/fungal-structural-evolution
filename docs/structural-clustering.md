# Exploratory structural similarity groups

The first clustering stage includes every model in the two frozen marker-model
snapshots: 13,153 AlphaFold and 5,121 ESMFold models, 18,274 total. Each receives
a unique alias linked to its original source, model ID/version, sequence hash,
coordinate path/hash, length and mean CA pLDDT. Identical sequences from different
prediction sources remain separate model records. Taxon and marker mappings
remain available in the pinned original snapshots; model clusters do not remove
species-specific variation from evolutionary inputs.

Foldseek runs its 3Di+AA alignment mode with sensitivity 7.5, up to 2,000
prefilter hits per query, E-value cutoff 1e-5, 80% aligned coverage in both
query and target and an approximate alignment-normalized TM-score threshold of
0.5. Greedy length-based single-step clustering and reassignment are requested.
The full command, installed executable hash/version and source receipts are
pinned in `metadata/frozen_structure_clustering_config.json`; option semantics
were checked against the installed executable's `easy-cluster -h` output.

pLDDT below 70 is masked for seeding only. This does not remove every uncertain
residue from alignment or enforce the six-residue/PAE gates used by paired
phylogenetic analyses. Clusters are exploratory similarity groups, not
orthogroups, validated homolog groups, functions or novel folds. Prefilter hit
limits may reduce recall. Representatives and members require explicit edge
validation; cluster membership does not imply all member pairs meet thresholds.
Threshold, prediction-source, domain-boundary and confidence sensitivities remain
required before interpretation. Full-chain coverage can miss valid domain-level
relationships, which require their own atlas analysis.

```bash
python scripts/cluster_frozen_marker_structures.py --output results/structural_clusters/frozen-marker-models-v1
```

The run is already active. It uses eight CPU threads, a 32 GB prefilter split
memory limit, 48 GB total memory planning and 100 GB output planning on the
existing host. The initial 1–48 hour forecast is broad pending production
timings; no paid infrastructure was provisioned. Inputs are symlinks to verified
source files, and temporary databases are retained for validation. The wrapper
rejects an existing output directory; inspect any partial failure before recovery
rather than starting a duplicate or overwriting results.

All 18,274 aliases and symlink targets passed readback, and the producer verified
every coordinate checksum. Foldseek database construction has started; there is
no completed clustering result yet. This is the full frozen **marker-model**
set, not the project's entire fungal proteome atlas.


## Initial clustering completed; membership verified

The run completed with 2,249 similarity groups across all 18,274 model records:
1,276 singleton groups and 302 groups containing both prediction sources.
Every Foldseek input database alias occurs exactly once. Independent readback
verified the complete output partition, representative self-membership, unique
annotated members, every retained provenance field, and agreement between raw
and annotated assignments. Source-specific taxon/marker links remain untouched.

These are algorithmic group counts, not numbers of orthogroups, novel folds or
independent evolutionary events. Alignment criteria have not yet been
independently checked for each representative/member pair. In particular,
confidence seeding does not certify all aligned residues. Edge validation,
threshold sensitivity, confidence/source analysis and domain-level clustering
remain pending. Large retained intermediate databases allow those follow-ups.

## All representative/member edges realigned

Recomputed all 16,025 non-self representative/member pairs in both directions
(32,050 directed alignments) without the clustering prefilter, retaining
backtraces and permissive output filters. Every requested alignment returned.
Scoring settings follow the clustering alignment log; this is the same Foldseek
implementation, not an independent geometry calculation.

Against the reported E-value, bilateral coverage and approximate alignment-TM
criteria, 15,043 pairs pass in both directions, 979 in neither, and three in only
one direction. The 982 non-bilateral pairs affect 200 original groups. Complete
identity, numeric-finiteness and threshold-decision readback passed; coverage
and TM-score failures dominate. The original memberships remain immutable and
cannot be treated as uniformly satisfying their representative/member criteria.

The audit also found 102 approximate alignment-normalized TM scores above one
(maximum 1.008). Those values are explicitly retained in the readback flags,
not silently clamped or presented as exact physical scores. Recalculation with
Foldseek `convertalis --exact-tmscore 1` is now running on all saved alignments.
Its command, input receipt/config hashes and resource estimate are saved in
`metadata/structure_cluster_exact_score_resource_config.json`. Recalculation
uses the same alignments, not a new search. Exact-score decision comparisons and
conservative reassignment remain pending.

```bash
python scripts/validate_structure_cluster_edges.py --clusters results/structural_clusters/frozen-marker-models-v1 --output results/structural_clusters/edge-validation-v1
```

The full realignment used eight CPU threads with 32 GB memory and 20 GB output
planning. Exact-score conversion reserves eight threads, 32 GB memory and 5 GB
output on the existing host, with a broad 0.1–12 hour forecast. No paid services
were used. Neither direction-specific failure nor singleton status establishes
biological dissimilarity, novelty or orthology.

## Marker and taxon provenance attached to original groups

All 2,249 original groups now have model/source/length/confidence summaries,
marker sets, taxon counts and major-lineage annotations. All 18,815
model–taxon–marker links exactly match the union of the frozen source tables;
per-group marker and taxon counts also passed independent readback. Twenty-three
groups contain multiple marker labels, and 95 contain outgroup taxa. These
patterns are review cues, not validated homology, functional or evolutionary
relationships. Shared domains, annotation differences and incorrect grouping
remain possible explanations for mixed marker labels.

The 200 groups with edges failing the approximate criteria retain explicit
counts and a pending-exact-review flag. Singleton groups have no edge test.
Every group keeps `orthology_status=not_established_by_structural_clustering`.
Original memberships have not been modified; exact-score conversion remains
running. A review script is prepared to compare all exact/approximate directed
pairs and reject changed alignment identities or out-of-range scores; it has
passed syntax compilation but has not yet run on completed exact-score output.

```bash
python scripts/annotate_structure_clusters.py --clusters results/structural_clusters/frozen-marker-models-v1 --edge-validation results/structural_clusters/edge-validation-v1 --output results/structural_clusters/provenance-annotations-v1
# Run only after exact conversion has a completion receipt:
python scripts/review_exact_cluster_scores.py --validation results/structural_clusters/edge-validation-v1 --exact results/structural_clusters/edge-exact-score-v1 --output results/structural_clusters/edge-exact-review-v1
```


### Exact-score review completed

Conversion with `--exact-tmscore 1` completed on the same 32,050 retained directed
alignments. All pair identities and fixed alignment fields match the preceding
review. Of 16,025 representative/member pairs, 14,769 pass both directions,
1,254 neither, and two one direction. There are 549 directed decision changes
and 275 pair classification changes relative to the approximate-score review.

Seventy-two directed alignment-normalized TM scores still exceed one (maximum
1.008); all other bounded score fields pass range checks. Direct table readback
preserves these values and marks their edges as failing review. The underlying
cause remains unresolved; the exact option alone does not resolve this issue.
No clamping, reassignment, orthology or novelty inference has been applied.
Original groups and their prior annotations remain unchanged and exploratory.

```bash
python scripts/review_exact_cluster_scores.py --validation results/structural_clusters/edge-validation-v1 --exact results/structural_clusters/edge-exact-score-v1 --output results/structural_clusters/edge-exact-review-v1
```

Receipts and the complete pair classification table are tracked under
`metadata/structure_cluster_exact_score_*`; full directed scores remain outside
Git under `results/structural_clusters/edge-exact-review-v1`.

### Conservative groups after edge review

`reviewed-groups-v2` retains 17,018 models in the 2,249 original representative
slots. Every retained nonself member (14,769) passes the fixed-alignment exact
score criteria in both directions. Representatives remain present without a
self-edge test. There are 1,315 slots with no retained nonself edge, compared
with 1,276 original singletons; this is a filtering result, not evidence for
39 newly discovered families.

The 1,256 deferred models have an empty derived group identifier: 1,220 fail
edge criteria and 36 have out-of-bounds scores. They are preserved in a separate
table with all original provenance. No alternative representative search or
novel-singleton interpretation is applied. All 18,815 taxon–marker links survive,
including deferred models, with explicit dispositions. Group summaries describe
retained membership only; original group annotations must not be silently
substituted for these summaries.

These are representative-centered similarity groups: passing links to the
representative do not guarantee similarity between every pair of members.
Orthology, domain-level robustness, confidence sensitivity and the cause of
out-of-range alignment-TM scores remain unresolved. The output covers the
frozen marker model collection, not the full proteome atlas.

```bash
python scripts/derive_reviewed_structure_groups.py --clusters results/structural_clusters/frozen-marker-models-v1 --review results/structural_clusters/edge-exact-review-v1 --annotations results/structural_clusters/provenance-annotations-v1 --output results/structural_clusters/reviewed-groups-v2
python scripts/audit_reviewed_structure_groups.py --groups results/structural_clusters/reviewed-groups-v2 --clusters results/structural_clusters/frozen-marker-models-v1 --exact results/structural_clusters/edge-exact-score-v1 --annotations results/structural_clusters/provenance-annotations-v1 --output results/structural_clusters/reviewed-groups-audit-v1
```

The readback independently recalculates eligibility from all 32,050 raw numeric
score rows and checks all 18,274 source models, membership partitions, group
model counts and the complete annotation multiset. It is not an independent
TM-score implementation. Version 2 adds a source-lineage consistency gate to
version 1 and produces identical four output tables; version 1 remains outside
Git. Receipts, summaries and deferred-model provenance are tracked in
`metadata/reviewed_structure_groups_*`; full membership and link tables remain
in the checksummed result directory.

### Coverage after conservative filtering

Coverage summaries include all 526 analysis-manifest taxa, including 67 with
no model in this frozen input collection. All 459 represented taxa retain at
least one model after filtering. This differs from paired-phylogenetic coverage,
which additionally requires eligible masks and enough taxa per marker.

| Source | Input models | Retained models | Deferred | Tested nonself members retained |
| --- | ---: | ---: | ---: | ---: |
| AlphaFold | 13,153 | 12,441 | 712 | 11,206 / 11,918 |
| ESMFold | 5,121 | 4,577 | 544 | 3,563 / 4,107 |

Overall retention is 94.6% and 89.4%, respectively. Representatives are retained
by definition and are reported separately; among tested nonself members the
fractions are 94.0% and 86.8%. These descriptive rates do not compare predictor
accuracy: source, protein length, confidence, taxon availability and representative
selection differ. Tables also stratify by marker, lineage, source, and source ×
length × mean predicted confidence. Mean confidence does not replace residue
confidence masks. Models are deduplicated within each summary; shared identical
models mean taxon/marker/lineage totals are not additive.

```bash
python scripts/summarize_reviewed_group_coverage.py --groups results/structural_clusters/reviewed-groups-v2 --manifest metadata/analysis_manifest.tsv --output results/structural_clusters/reviewed-coverage-v1
```

Independent pandas count readback checked all 677 taxon, marker, lineage and
source summary rows and the complete 526-taxon universe. The additional
source-length-confidence table was not included in that independent readback.
Tables and receipts are tracked in `metadata/reviewed_structure_coverage_*`.
Lineage rows describe available models; the taxon table preserves absent inputs.

### Source and alignment traces explain the possible score overshoot

Retrieved and hashed scoring source at the exact reported binary revision
`e3fadcd07f971e864c094ac4f3a78bf4ed845e07`. The
[alignment-score conversion](https://github.com/steineggerlab/foldseek/blob/e3fadcd07f971e864c094ac4f3a78bf4ed845e07/src/strucclustutils/structureconvertalis.cpp#L996)
passes `min(qEnd-qStart, tEnd-tStart)` as normalization length. All 32,050
exported CIGARs reconstruct inclusive endpoints, so that expression omits one
position from the shorter span. The
[exact-score implementation](https://github.com/steineggerlab/foldseek/blob/e3fadcd07f971e864c094ac4f3a78bf4ed845e07/src/commons/TMaligner.cpp#L107)
uses normalization in the final distance scale and score optimization.

The full trace audit finds 7,398 alignments where the implemented denominator
is smaller than the number of matched positions. All 72 observed scores above
one are in this set and satisfy the corresponding upper bound
`matched_positions / implemented_denominator` (allowing printed precision).
For example, M000926→M011904 is `122M` over positions 1–122, denominator 121,
with score 1.008. This supports an off-by-one normalization explanation, but
has not yet been confirmed by running a corrected binary. Merely multiplying
scores by a length ratio would not repeat the altered optimization. Deferred
members remain deferred pending corrected rescoring and threshold sensitivity.

```bash
python scripts/audit_cluster_score_normalization.py --traces results/structural_clusters/edge-backtrace-review-v1 --review results/structural_clusters/edge-exact-review-v1 --source data/software_source/foldseek-score-review-v1 --output results/structural_clusters/normalization-audit-v1
```

The trace export command, source URLs/hashes and complete 72-case table are
tracked in `metadata/cluster_score_normalization_*`. Source and full trace
tables are outside Git. No installed executable or original result was altered.
