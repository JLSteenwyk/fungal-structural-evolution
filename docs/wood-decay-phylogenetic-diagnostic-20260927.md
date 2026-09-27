# Placement uncertainty for reviewed wood-decay classifications

All 1,049 edges of each of two completed species trees were evaluated under
five explicit coding scenarios. The primary scenario uses four reviewed
white-rot labels (0), two brown-rot labels (1), and leaves both uncertain taxa
and all other taxa unknown (520 unknown tips). Four sensitivities assign each
possible binary combination to Jaapia and Botryobasidium. These assignments
explore a model's sensitivity; they do not change the curated source labels or
establish that a binary decay model is biologically adequate.

Both trees give the following conditional results:

| Jaapia / Botryobasidium coding | Minimum changes | Required edges | Optional edges |
| --- | ---: | ---: | ---: |
| Unknown / unknown | 2 | 0 | 10 |
| White / white | 2 | 0 | 5 |
| White / brown | 2 | 0 | 6 |
| Brown / white | 2 | 1 | 3 |
| Brown / brown | 2 | 1 | 4 |

A required edge changes under every optimal assignment. Optional edges change
under some but not all optima; their assignments are jointly constrained.
The primary scenario therefore does not localize either minimum change to a
unique edge. A required edge appears only when the uncertain Jaapia label is
forced to brown. That split contains F104355, F202697, F38799 and F5364 on one
side in both trees. It is a coding-dependent result, not a validated origin.
Two brown-rot tips and a minimum of two undirected changes do not establish two
independent gains or two usable structural contrasts. Rooting, tree-ensemble
uncertainty, sparse observations and the trait model remain unresolved.

The inside/outside recurrence computed all four conditional endpoint costs for
every edge. Independent integer network maximum flow verified **41,960 costs**
across all **10,490 edge/scenario/tree rows**, as well as identities, optimal
assignments and summaries. No producer recurrence is imported by the checker.
The existing exhaustive small-tree tests support the shared recurrence; the
new coding and complete production outputs were checked by the separate flow
algorithm. Both executions completed successfully.

Reproduction (use fresh output paths):

```bash
python scripts/map_wood_decay_optimal_edges.py --plan metadata/wood_decay_optimal_edges_plan_20260927.json
python scripts/readback_wood_decay_optimal_edges.py --plan metadata/wood_decay_optimal_edges_plan_20260927.json --output metadata/wood_decay_optimal_edges_readback_20260927.json
```

Plans pin the reviewed evidence and source trees. Full edge tables are outside
Git at `results/ecology/wood-decay-optimal-edge-states-20260927-v1`.
The source-bound completion record is
`metadata/wood_decay_optimal_edges_completed_20260927.json`, and the ten-row
summary is `metadata/wood_decay_optimal_edges_summary_20260927.tsv`.
Planning allowance: one CPU, 4 GiB RAM, 1–10 minutes per stage, no GPU.
The [coverage and source review](ecology-source-review-20260927.md) retains the
predictor imbalance and uncertain classifications needed to interpret this result.

## Complete bootstrap mapping launched

The same five coding scenarios are now running across all 2,000 saved bootstrap
trees, yielding 10,000 tree/scenario combinations. Every edge's four conditional
endpoint costs are retained in hashed compressed per-tree files. Per-split
summaries distinguish presence from required, optional and excluded changes;
absence from a tree is not treated as no change. The new scripts leave the
existing mycorrhizal-state bootstrap outputs intact.

Before queuing the full readback, the independent network-flow kernel checked
all 20,980 endpoint costs of one complete native bootstrap tree across the five
scenarios in 17.16 seconds. This validates the new scenario wiring, not the full
ensemble. The queued checker requires terminal producer success and then
recomputes every cost across all 2,000 trees; full verification is pending.
Expected full scope is 10,490,000 edge rows and 41,960,000 constrained costs.

Producer: `scripts/map_wood_decay_bootstrap_edges.py`, plan/launch
`metadata/wood_decay_bootstrap_{plan,launch}_20260927.json`.
Readback: `scripts/readback_wood_decay_bootstrap_edges.py` with independent kernel
`scripts/wood_decay_bootstrap_flow_readback.py`, plan/launch
`metadata/wood_decay_bootstrap_readback_{plan,launch}_20260927.json`.
Native check: `metadata/wood_decay_bootstrap_native_checks_20260927.json`.
Outputs: `results/ecology/wood-decay-bootstrap-{edges,readback}-20260927-v1`.

The producer has one CPU, 8 GiB RAM, no swap, a 5 GiB output allowance and a
5–120 minute planning range. Readback has four CPUs, 8 GiB RAM, no swap, a 5 GiB
output allowance and a 1–24 hour planning range after production. These are
resource allowances, not completion forecasts. GPU inference remains paused.
Bootstrap frequencies remain conditional topology sensitivity, not posterior
transition probabilities, independent origin counts or structural-effect tests.

## Bootstrap uncertainty summaries verified

The full independent network-flow audit completed all 2,000 trees and five
codings: 10,490,000 edge rows and 41,960,000 endpoint-constrained costs. All 2,000
per-tree proof hashes were checked. Summary production and independent full
readback also completed: 12,515 split rows, 74 metric-distribution rows and
8,000 paired coding-sensitivity rows. All four services terminated successfully.

For the primary coding, both 1,000-tree ensembles have minimum change count two
on every tree. However, every tree has **zero edges required to change in all
optimal mappings**. This does not mean no changes occurred: the minimum changes
can be assigned to different edges in alternative optima. Primary optional-edge
counts range from nine to eleven.

Only the scenarios forcing Jaapia to brown produce a required-change edge. In
each ensemble, 996 of 1,000 trees require the split whose canonical side is
`F104355;F202697;F38799;F5364`; that split is present in all 1,000 trees. The four
remaining trees have no required edge. Assigning Botryobasidium either white or
brown does not alter that 996/1,000 result. The 99.6% frequency is conditional on
an uncertain trait assignment; it does not establish a brown-rot origin or an
independent ecological transition. The primary labels remain unchanged.

This diagnostic therefore does not supply confidently localized transitions for
the primary wood-decay/structural-change association test. Sparse trait coverage,
alternative optimal mappings and uncertain assignments remain substantive
limitations. The broader ecological aim remains open.

`summarize_wood_decay_bootstrap_uncertainty.py` preserves split-presence and
change frequencies with explicit all-tree and present-only denominators. Paired
sensitivity ranges compare alternative codings on the same tree and are not
confidence intervals. `readback_wood_decay_bootstrap_summary.py` independently
reconstructed every summary field with full joins and grouped counts.

Output: `results/ecology/wood-decay-bootstrap-summary-20260927-v1`.
Completion: `metadata/wood_decay_bootstrap_completed_20260927.json`.
Independent summary readback:
`metadata/wood_decay_bootstrap_summary_readback_20260927.json`.
Reproduce the completion binding with
`scripts/record_wood_decay_bootstrap_completion.py` after the recorded full jobs
have terminated successfully. This completion covers the diagnostic computation,
not a test of structural effects or completion of the ecological project aim.
