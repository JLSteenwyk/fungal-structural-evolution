# Coalescent and concatenated candidate-tree comparisons

The full comparison is queued after the independent numerical closure of all
30 ASTRAL candidates. It compares candidate relationships across gene-tree
alignment/support settings and against every original or exactly matched
taxon projection of the four concatenated PMSF runs. These comparisons inform
the tree uncertainty needed for subsequent structural-rate and reconciliation
analyses. They do not select an accepted species tree or biological root.

## Complete design

Each of five cohorts has six coalescent candidates (two gene alignments by
three support settings) and eight concatenated reference views (four original
fits by ML/consensus). Every coalescent candidate is compared against every
same-cohort reference: 240 comparisons. All pairs among the six coalescent
candidates contribute another 75 comparisons. The total is 315 comparisons
across 70 views; no candidate, reference or conflicting split is selected out.

The primary cohort retains 526 entries, including 25 outgroups. The four closed
retained cohorts contain 525/522/524/503 entries and 25/23/25/25 outgroups.
Reference projections retain the original fitted branch-path sums and
recomputed bootstrap frequencies; they are not new subset inference.
The 16 actual supported subset PMSF refits and their comparisons remain
separate requirements.

For every cohort, the output includes every view's presence or absence at
every split in the union, all matching and unique split counts, Robinson–Foulds
distance, every incompatible unique-split pair and a canonical four-taxon
witness. Pairwise distances use identical taxon sets. All 70 declared
ingroup/outgroup boundaries are checked as unrooted splits and retain their
actual membership counts; their presence does not establish biological rooting.

## Support and branch units

Separate descriptive support screens are preserved:

- Coalescent candidates: native local posterior at least 0.95, with effective
  gene count retained separately.
- Original concatenated ML: SH-aLRT at least 80 and empirical UFB at least 95%.
- Original consensus: empirical UFB at least 95%; SH-aLRT unavailable.
- Projected references: recomputed empirical UFB at least 95%; SH-aLRT unavailable.

These measures are not interchangeable probabilities or calibrated global
confidence. A pair passing both applicable screens is a descriptive overlapping
conflict, not an independent evolutionary event or evidence of a particular
cause of discordance. Unavailable support remains unavailable.

Lengths are retained with explicit units: MAP coalescent units, amino-acid
substitutions per site, or sums of original amino-acid branch lengths on
projected paths. No numerical length comparison is made between methods.
Likelihoods from different site universes are not compared.

## Verification and completed first-candidate evidence

The producer uses integer split masks and four nonempty intersection cells to
identify incompatible splits. A separate reader parses every actual native
tree in DendroPy, prunes original reference trees to the exact retained cohort,
reconstructs every split and branch length, and checks RF using DendroPy's
tree comparison and incompatibility using its bipartition implementation
([Sukumaran and Holder 2010](https://doi.org/10.1093/bioinformatics/btq228);
[Moreno et al. 2024](https://doi.org/10.21105/joss.06943)).
It checks every exported presence cell, metric/unit, support disposition,
boundary and quartet witness. Full source/output hashes and both original
producer/reader completion journals are required for closure.

The full five-cohort/70-view/315-pair synthetic software grid passed, including
exact/below/above local PP 0.95, UFB95 and SH80 boundaries. The independent
reader rejected 13 changed support, metric, topology-count, witness or boundary
exports. These are software contracts, not biological results or a pilot.

The actual first full 526-taxon profile/uncontracted candidate was compared
with all eight original reference views. All eight comparisons, 5,814 presence
cells, 1,988 overlapping incompatible split pairs and nine role boundaries
passed independent raw-tree reconstruction. Closure verified 1,475 bindings
and both original completion journals. [Completed named-case locator](../metadata/coalescent_comparison_preflight_completed_20261002.json).

This candidate shares 431–438 of 523 internal splits with each reference;
the RF distances are 170–184. Across the candidate and all eight references,
417 splits are shared. Each pair has 42–70 incompatible pairs passing the
two applicable support screens. These are dependent comparisons conditional
on inferred genes and candidate trees; they do not establish biological
discordance, model adequacy or a preferred estimate. The primary role boundary
is present in this candidate (local posterior 0.9982808647585494, effective
genes 107); this remains an unrooted diagnostic. The other 29 candidates are
not certified by this named-case comparison.

## Execution and remaining work

Resources were estimated before launch: two CPU/16 GiB/no swap/one BLAS
thread for each serial producer, reader and closure; 4 GiB output allowance
and 100 GiB disk reserve. The actual first eight comparisons took 8.506 CPU
seconds for the producer and 9.403 for the reader. The full conservative
candidate-check bound is 315 times 523 squared, before removing shared splits.
The 0.02–2 hour per-stage planning range excludes dependency waits and is not
a guaranteed ETA. No inference is restarted, and no GPU or paid resource is used.

```bash
python scripts/compare_coalescent_species_trees.py \
  --plan metadata/coalescent_tree_comparison_plan_20261002.json
python scripts/readback_coalescent_tree_comparisons.py \
  --plan metadata/coalescent_tree_comparison_plan_20261002.json \
  --output results/phylogeny/full-coalescent-reference-tree-comparisons-20261002-v1/readback.json
python scripts/close_coalescent_tree_comparisons.py \
  --plan metadata/coalescent_tree_comparison_completion_plan_20261002.json
```

Commands identify the frozen queued runs; do not start duplicate runs.
Reproduction uses new plans/output locations. Full tables/conflict records,
trees and hash dictionaries stay outside Git. The future small locator is
`metadata/coalescent_tree_comparisons_completed_20261002.json`.

The [current runtime checkpoint](../metadata/project_runtime_checkpoint_20261002_v3.json)
records 17/30 native coalescent cases complete, 25 live pipeline handles and
six original scientific/retrieval jobs; new structural comparisons have
processed 162,816/298,848 states at its observation time. Full coalescent
numerical readback/closure and these 315 comparisons remain queued. Gene-tree
estimation uncertainty, marker dependence/selection, taxon identity, MSC and
sequence-model adequacy, roots, accepted dating and reconciliation remain open.
All eight evolutionary aims remain incomplete. GPU prediction remains paused.
