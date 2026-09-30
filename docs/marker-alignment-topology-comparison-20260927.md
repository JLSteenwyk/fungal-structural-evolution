# Marker alignment topology comparison

**September 30 status:** the full 125-marker MAFFT inference, source-input/support
audit, common-taxon topology comparison and independent readback have completed.
The earlier preparation and queued-status text below is retained as history;
completed results and the follow-up support analysis are recorded at the end.

The full 125-marker MAFFT inference and input/support audit must complete before
comparison with the audited profile-alignment trees. Thirty markers differ in
retained taxon membership (149 additional MAFFT marker–taxon entries); the
comparison therefore verifies each emitted tree's exact coverage-derived tips
and restricts its splits to the shared tips. Duplicate projected splits are
collapsed and trivial splits omitted. This is equivalent to pruning trees, not
refitting them; it retains effects that different original taxon sampling had
on inference. Zero-length internal edges remain represented.

`scripts/compare_marker_alignment_topologies.py` requires both completed full
support snapshots and their independent audit receipts. It validates tree hashes,
case universes and membership before writing all shared/lost/gained split rows
and per-marker unrooted Robinson–Foulds distances. No likelihood or support test
is implied by topology differences. Independent readback of the eventual full
cross-method output remains required.

Verification completed using all 125 existing profile trees:

- Full identity comparison: zero differences, 59,315 internal split rows, all
  matching the independently audited support table.
- Independent DendroPy pruning versus Bio.Phylo split restriction: all 125 trees
  under both full-tip and all-outgroup-excluded conditions, 250 comparisons and
  116,039 internal split checks, with exact agreement.

Evidence is archived in `metadata/marker_topology_identity_readback_20260927.json`
and `metadata/marker_topology_projection_readback_20260927.json`; the latter is
reproduced with `scripts/check_marker_topology_projection.py`. These checks
validate the comparison machinery on existing data; they do not establish
results for unfinished MAFFT trees.

After the complete MAFFT audit succeeds, run:

```bash
python scripts/compare_marker_alignment_topologies.py \
  --first-trees results/phylogeny/marker-gene-trees-v2 \
  --first-support results/phylogeny/marker-tree-support-complete-v1 \
  --first-audit results/phylogeny/marker-tree-support-complete-readback-v1 \
  --first-coverage profile \
  --second-trees results/phylogeny/mafft-marker-gene-trees-20260927-v1 \
  --second-support results/phylogeny/mafft-marker-support-complete-20260927-v1 \
  --second-audit results/phylogeny/mafft-marker-support-readback-20260927-v1 \
  --second-coverage mafft \
  --coverage results/phylogeny/marker-alignment-coverage-20260927-v1 \
  --output results/phylogeny/marker-alignment-topology-comparison-20260927-v1
```

The cross-method comparison and independent readback are now queued after the
complete MAFFT input/support audit; neither cross-method result has completed.


The independent readback script is
`scripts/readback_marker_alignment_topologies.py`. It accepts the same source
arguments plus `--comparison` and a new JSON `--output`. All 125 marker summaries
and 59,315 split rows passed its full DendroPy-based identity-batch readback.
Evidence: `metadata/marker_topology_identity_independent_readback_20260927.json`.
The earlier 250-condition pruning check separately covers taxon removal.

Automatic handoff uses `scripts/advance_marker_alignment_topology_comparison.py`
and `metadata/marker_alignment_topology_handoff_plan_20260927.json`. The waiting
service is confirmed live with one CPU, 8 GiB and no swap. It pins the scripts
and completed sources, requires the entire MAFFT audit receipt, then compares
and independently verifies all markers. Allowance after dependencies is 0.1–2
hours and 1 GiB output. Comparison results still await those dependencies.

## Full completed results, September 30

All 125 MAFFT marker trees completed native inference with return code zero.
The independent input/support readback checked 30,575,628 retained alignment
characters and 59,464 internal splits. The complete cross-method readback then
verified 80,512 shared/lost/gained split rows and every marker summary using
DendroPy pruning separately from the producer's split restriction. The
[completion record](../metadata/marker_alignment_topology_completed_20260930.json)
rechecks all tree/input and receipt hashes and the captured process journals;
the three resumed services have finished. Earlier interruption records remain
preserved.

All 125 marker topologies differ after projection to shared taxa. Thirty
markers had different original retained taxon memberships; the remaining 95
had the same tips. Across markers, the unrooted Robinson–Foulds distance ranges
from 92 to 700, with median 342. Dividing by twice the shared taxon count minus
six gives a median normalized distance of 0.3540 (range 0.0979–0.7261).
The median Jaccard similarity between internal-split sets is 0.4771
(range 0.1587–0.8217). These summaries retain zero-length edges and do not
condition on support. They do not identify significant conflicting branches,
a preferred alignment method or a biological cause of discordance. Pruning
does not remove the effect that original taxon sampling had on inference.

The full output tables remain outside Git at
`results/phylogeny/marker-alignment-topology-comparison-20260927-v1`;
their checksums are bound through the completion record. Reproduction uses the
unchanged commands above and a new output path, followed by
`scripts/readback_marker_alignment_topologies.py`.

## Full supported-conflict sensitivity completed

The new [full-scope plan](../metadata/alignment_supported_conflict_plan_20260930.json)
uses both audited homogeneous full-taxon guides and all 125 markers under both
alignment methods. For every guide internal edge it retains exact matching
SH-aLRT support, the strongest incompatible split and its witness, and explicit
concordant/conflicting/unresolved/uninformative classifications at cutoffs 80
and 95. The original profile diagnostic is reused without alteration; the full
MAFFT diagnostic is new. A source-linked paired classification table is produced
only after both independent checks pass.

The new checker independently searches all 130,750 marker/guide-edge conflict
maxima under each method using Python sets; the producer uses bit-vector
intersections. It verifies all 261,500 cutoff rows per method, including equal
stored maxima across both cutoffs, exact supports, witnesses, complete edge
grids and all 2,092 summaries. The older checker searched maxima on only five
edges per marker/guide, so its proof remains distinct. The
[validation record](../metadata/supported_conflict_exhaustive_fixture_checks_20260930.json)
covers all four support/coverage states and rejection of rehashed false maxima,
including an inconsistent value under just one cutoff. This is machinery
validation, not a separate biological pilot or proof of full-source completion.

The [live launch](../metadata/alignment_supported_conflict_launch_20260930.json)
uses one CPU, 16 GiB memory, zero swap and no GPU or paid infrastructure.
The prelaunch estimate allows 2 GiB output and 0.5–24 active hours for full
diagnostic production, two exhaustive independent checks and the serialized
join. This broad planning allowance was not a measured ETA or timeout. The
service finished at 10:58 EDT on September 30, using 189.554 CPU seconds,
363.5 MiB peak memory and zero swap. All full stage receipts and their source
hashes passed the [completion check](../metadata/alignment_supported_conflict_completed_20260930.json).

Both exhaustive checks passed all 261,500 cutoff rows and all 130,750 independent
conflict maxima per method. A further independent pandas one-to-one outer merge
reconstructed every paired classification and all 52 transition-count rows from
the audited sources. This complete
[join readback](../metadata/alignment_supported_conflict_join_completed_readback_20260930.json)
also produces 2,092 branch-level sensitivity rows and eight summaries stratified
by original retained tip membership. The 95 same-tip markers and 30 different-tip
markers remain separate; differences cannot be attributed solely to alignment.

| Full guide | SH-aLRT cutoff | Marker–edge cells | Changed classification | Concordance ↔ conflict reversal | Concordant under both | Conflicting under both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Profile | 80 | 65,375 | 13,621 | 2,675 | 23,875 | 16,484 |
| MAFFT | 80 | 65,375 | 13,719 | 2,690 | 24,073 | 16,401 |
| Profile | 95 | 65,375 | 11,735 | 404 | 16,376 | 6,132 |
| MAFFT | 95 | 65,375 | 11,790 | 402 | 16,520 | 5,983 |

Across the full 261,500-cell grid, 50,865 classifications change. Changed
classifications also include transitions to/from unresolved or uninformative
coverage, and therefore are broader than direct concordance/conflict reversals.
The guides and cutoffs overlap, and cells can share markers, taxa and collapsed
paths; the counts must not be pooled as independent evolutionary events.

The branch and membership tables remain at
`results/model_validation/alignment-supported-conflict-readback-20260930-v1`,
bound by the versioned completion record. Reproduce the full paired readback
with `python scripts/readback_alignment_supported_conflict_sensitivity.py
--plan metadata/alignment_supported_conflict_join_readback_plan_20260930.json` using
a new output path and corresponding new pinned plan. Completed outputs are
preserved; rerunning native marker inference is unnecessary for these tables.

These support-cutoff summaries are descriptive. SH-aLRT is neither bootstrap
probability nor a gene concordance factor. Marker coverage can collapse guide
paths, and rows are not independent evolutionary events. Species-tree model
adequacy, supported mixture inference, root/dating sensitivity and biological
interpretation remain required before structural changes can be assigned to
accepted evolutionary branches.
