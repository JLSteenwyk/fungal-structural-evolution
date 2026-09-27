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

## Bootstrap uncertainty summaries queued

The full 2,000-tree, five-coding network-flow audit remains active. After its
successful terminal state, `summarize_wood_decay_bootstrap_uncertainty.py` will
produce split-level presence and change frequencies with explicit all-tree and
present-only denominators, plus tree-level metric distributions. It also computes
8,000 paired tree/metric rows contrasting the primary unknown coding with the
minimum and maximum over the four uncertain-taxon assignments on the same tree.
These are coding sensitivity bounds, not confidence intervals.

`readback_wood_decay_bootstrap_summary.py` waits for successful summary completion
and independently reconstructs every field using pandas outer joins and grouped
counts. Both stages use one CPU, 4 GiB, no swap, with 1–15 minutes estimated per
stage after its prerequisite. Launch identities and source/script hashes are in
`metadata/wood_decay_bootstrap_summary_launch_20260927.json` and
`metadata/wood_decay_bootstrap_summary_readback_launch_20260927.json`.

No bootstrap summary is yet accepted. The frozen primary trait assignments remain
unchanged; changes in optimal mappings do not establish rooted origins, independent
replication, transition probabilities, or effects on protein structure.
