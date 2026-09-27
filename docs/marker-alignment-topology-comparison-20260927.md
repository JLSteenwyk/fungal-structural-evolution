# Marker alignment topology comparison

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
