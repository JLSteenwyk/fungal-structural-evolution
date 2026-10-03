# Structural observations across candidate phylogenies

This audit maps the complete qualified AlphaFold and ESMFold marker memberships
onto **all 70 closed candidate tree views**, covering five declared taxon cohorts.
It keeps all 125 marker slots and every original internal branch. Its purpose is
to show where missing structural observations prevent a branch from remaining
distinct when the tree is pruned. It does not estimate evolutionary rates.

The primary panel contains 526 working taxon entries: 501 fungal entries and
25 outgroups. Species boundaries remain under review. Sensitivity cohorts retain
their declared membership; two have 23 outgroups. These exclusions do not change
the primary panel. Each cohort has six coalescent and eight concatenated views.
Projected concatenated sensitivity views are prunings of existing fits, not
newly inferred subset phylogenies. Native subset refits remain separate work.

## Inputs and identity

The full sources are `paired-inputs-afdb-recovered-20260925-v1` and
`paired-inputs-esmfold-all-completed-20260922-v1` under `results/phylogeny`.
Their original complete readbacks qualify 31,134 AlphaFold marker–taxon cells
across 125 ready markers and 339 taxa, and 23,799 ESMFold cells across 122 ready
markers and 294 taxa. All 125 ESMFold slots remain in this audit, including
eligible observations in the three families with fewer than four eligible taxa.
Sources remain separate; their complete eligibility and paired-alignment masks
are verified before projection. Existing prediction-source calibration is not
established by these checks.

Each usable taxon requires at least `max(50, ceil(0.3 * original_marker_columns))`
observed paired columns. Existing qualification requires canonical mapped amino
acids, valid native structural states, six-feature pLDDT at least 70 and
directional PAE at most 10. The audit does not relax those requirements.

The table retains frozen names and supplies display/provenance names separately.
It applies the exact assembly-bound review of O2700040 as **Branchiostoma
belcheri (hybrid-derived principal haplotype)**. Hybrid specimen origin does not
establish mixed-parental assembly ancestry. No original manifest, query identifier,
tree fit or native input is changed. Manifest lineage bins are descriptive
labels, not established monophyletic groups or biological species counts.

## Complete projection and independent readback

There are 17,500 view × predictor × marker cases, 36,190 original view-internal
branches and 9,047,500 branch-projection entries. Each original branch is assigned
one of five states:

| Code | State | Meaning |
|---|---|---|
| 0 | `insufficient_family_taxa` | Fewer than four observed taxa remain in this view/family/source |
| 1 | `no_observations_on_one_side` | No observed taxon remains on one side of the original split |
| 2 | `terminal_projection` | Exactly one observed taxon remains on one side |
| 3 | `unique_internal_projection` | At least two observed taxa remain on each side and only this original internal branch contributes to the projected split |
| 4 | `shared_internal_projection` | Multiple original internal branches contribute to the same projected split |

Unrooted bipartition orientation uses tip count followed by stable numeric mask
order, solely for bookkeeping. It does not accept a root. For each retained
internal split, the producer sums the lengths of its contributing original
internal branches. These remain in each source view's units: amino-acid
substitutions per site or native MAP coalescent units. They are not physical
displacement, time, newly fitted structural change or comparable rate scales.
Terminal, lost-side and insufficient-family entries retain **NaN**, not zero,
for the unclaimed path length. Observed tip counts remain explicit.

The independent reader parses each raw Newick tree with DendroPy, prunes it
first to its declared cohort and then to each exact source/marker membership,
and compares every saved entry against the actual graph's splits and path
lengths. It checks all five array names, dtypes, shapes, complete marker/taxon
axes, summary rows and original branch identities. Shared input provenance and
source qualification are reused; the pruning/path oracle is independent of the
producer's bitmask grouping algorithm. All source and artifact hashes are checked
before and after readback.

Every stored `view_XX.npz` has dimensions `(2 predictors, 125 markers,
original_internal_branches)`. Fields are `status` (uint8), `observed_side`,
`observed_complement`, `path_internal_branch_count` (uint16) and `path_length_sum`
(float64). Original branch masks, units and view membership live in
`input_axes.json`; the 526-row identity/source table and complete 17,500-row
summary remain alongside the arrays.

## Execution evidence and current status

The software gate passed every one of the 512 retained subsets of a synthetic
nine-tip unrooted tree against actual graph pruning. All five classes occurred;
11 altered value/serialization cases were rejected. All real source grids were
checked. This is software qualification, not a fungal pilot.

The full producer finished with exit zero and reopened every serialized array.
**Independent readback is running; full computational closure is pending.**
[Original readback checkpoint](../metadata/structural_marker_tree_projection_execution_checkpoint_20261003_v1.json).
The separate descriptive summary and figure script is prepared but must wait
for the complete readback receipt. It includes all 70 views, both sources and
every declared cohort; repeated views and marker entries are not independent
statistical observations. A distinct projected branch is not proof of statistical
estimability, sufficient information or an accepted sequence/structure model.

The prelaunch resource plan allows two CPUs, 8 GiB memory, no cgroup swap and
one BLAS thread. Each stage has an 8 GiB address-space cap, 1,800 CPU seconds,
2,400 wall seconds and a 512 MiB per-file limit. Planning allows 1 GiB output
and requires 64 GiB free disk and 16 GiB available RAM. These are bounds and an
uncalibrated planning range, not an ETA. Original software and producer stages
used 17.31 and 44.90 child CPU seconds, with measured child peak RSS about
264 and 254 MiB. Raw manager memory counters are preserved separately and are
not used as a native/group memory requirement. No GPU prediction, native sampler,
paid infrastructure or existing-job restart was launched by this audit.

## Reproduction

Use the repository Python environment with NumPy and DendroPy. Run each command
inside the documented resource limits, using new output paths. Launched sources,
attempts and receipts must remain unchanged.

```bash
python scripts/check_structural_marker_tree_projection.py \
  --output NEW_SOFTWARE_DIRECTORY --receipt NEW_SOFTWARE_RECEIPT
python scripts/project_structural_markers_on_candidate_trees.py \
  --output NEW_PROJECTION_DIRECTORY --receipt NEW_PRODUCER_RECEIPT \
  --resources metadata/structural_marker_tree_projection_resources_20261003_v1.json \
  --software NEW_SOFTWARE_RECEIPT
python scripts/readback_structural_marker_tree_projection.py \
  --producer NEW_PRODUCER_RECEIPT --output NEW_READBACK_RECEIPT
python scripts/summarize_structural_marker_tree_projection.py \
  --producer NEW_PRODUCER_RECEIPT --readback NEW_READBACK_RECEIPT \
  --output NEW_SUMMARY_DIRECTORY --receipt NEW_SUMMARY_RECEIPT
```

`run_structural_projection_stage.py` records original PID/create times, child
exits, logs, actual cgroup limits and resource use for the first three stages.
`close_structural_marker_tree_projection.py` additionally checks full hashes,
original wait-tool completion evidence and matching invocation startup/resource
journals. It never treats collected systemd default-success properties as
completion evidence. Full accepted species/reconciled trees, branch-specific
sequence/structural estimates, calibration and all eight evolutionary aims
remain incomplete.
