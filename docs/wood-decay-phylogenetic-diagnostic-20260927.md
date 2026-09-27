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
