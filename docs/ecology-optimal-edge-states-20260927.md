# Locations allowed by minimum-change ecological reconstructions

Mapped every edge on both completed profile-alignment PMSF trees under the
previously frozen 32-species binary classification. The 24 ectomycorrhizal and
seven asymbiotic/saprotrophic classifications retain their earlier conditional
coding; orchid symbiosis and the other 495 tips are unconstrained. The thirteen
newly curated species are not newly coded as ECM absence. A second scenario
leaves Ramaria rubella unknown because of its documented source uncertainty.

For each edge, an inside/outside dynamic program calculates the global minimum
score under each of its four endpoint-state combinations. It identifies edges
that change in every optimum, in some optima, or in none. This evaluates all
optimal assignments without choosing one arbitrary ancestral reconstruction.

| Conditional coding | Minimum total changes | Required-change edges | Optional-change edges | Edges unchanged in every optimum |
|---|---:|---:|---:|---:|
| Original 32-species coding | 7 | 3 | 19 | 1,027 |
| Ramaria unknown | 6 | 2 | 19 | 1,028 |

Both tree alternatives give these same counts. Required-change edges in the
first scenario are terminal edges for Botryobasidium botryosum (F264124),
Sphaerobolus stellatus (F68786) and Ramaria rubella (F113071). Leaving Ramaria
unknown removes its required change. These are **conditional constraints on
an undirected equal-cost reconstruction**, not experimentally established
transitions, three independent origins, or inferred biological gain/loss
polarities. The stored parent/child direction follows Newick traversal and
must not be interpreted as a biologically justified ancestral direction.

Optional edges are not independent choices: assignments across edges must
jointly attain the global minimum. In particular, the 19 optional edges do not
represent 19 events. Sparse and uncertain trait observations can constrain a
minimum change while leaving its location or direction unresolved. Identical
counts on the two trees do not establish support across bootstrap trees,
alignment alternatives, trait-coding schemes or evolutionary models.

## Full verification and reproducibility

Exhaustive enumeration passed 243 small-tree state patterns, including unknown
tips, alternative resolutions and a four-way polytomy. An independent integer
network-flow algorithm then checked all 16,784 endpoint-constrained costs
across 4,196 production edge rows, plus split identities, optimal assignment
sets, unknown-tip counts and every summary. No dynamic-programming helper was
imported by the production readback. Infeasible endpoint assignments have blank
costs, rather than interpreting the implementation's finite penalty as a cost.

```bash
python scripts/check_ecology_optimal_edge_states.py
python scripts/map_ecology_optimal_edges.py --plan metadata/ecology_optimal_edges_plan_20260927.json
python scripts/readback_ecology_optimal_edges.py --plan metadata/ecology_optimal_edges_plan_20260927.json --output metadata/ecology_optimal_edges_readback_20260927.json
```

Use fresh output paths for reruns. Full results are outside Git under
`results/ecology/optimal-edge-states-20260927-v1`. The plan, producer receipt,
full network-flow proof, script hashes, versions and reproduction commands are
in `metadata/ecology_optimal_edges_*_20260927.json`. The computation uses local
CPU only. The original tree files and their independent proofs are pinned.

Further work must evaluate edge-location uncertainty across the relevant tree
ensemble, review biological rooting and explicit trait models, and establish
usable replicated transitions with matched structural coverage. This diagnostic
advances phylogenetic integration but does not complete ecological testing.

## Full bootstrap edge mapping running

Started the same optimal-edge calculation across all 2,000 saved ultrafast
bootstrap trees (1,000 per guide) under both frozen coding scenarios. All 526
tip identities are checked per tree. Compressed per-tree checkpoints retain
every split, endpoint identity, status and four conditional costs. Summaries
separate the number of trees containing a split from counts of required,
optional and excluded changes, so an absent split is not treated as an observed
unchanged edge. The canonical split side is root-independent; state ordering
remains a traversal convention.

The pre-run allowance is one CPU, 4 GiB RAM, no swap, 2 GiB output and 2–30
minutes. The launched systemd unit enforces the CPU and memory allowances.
PID, creation time, exact command and plan hash are recorded in
`metadata/ecology_bootstrap_edges_launch_20260927.json`; the frozen plan is
`metadata/ecology_bootstrap_edges_plan_20260927.json`. Outputs are under
`results/ecology/bootstrap-optimal-edge-states-20260927-v1`. Resume reuses only
shards with matching plan/source hashes and verified compressed artifact hashes.

```bash
python scripts/map_ecology_bootstrap_edges.py --plan metadata/ecology_bootstrap_edges_plan_20260927.json
```

The validated recurrence is unchanged. Production completion and independent
bootstrap verification remain pending. These frequencies assess sampled tree
topology sensitivity conditional on one sequence alignment and fixed trait
coding; they are not posterior transition probabilities or independent origins.

## Full bootstrap network-flow verification queued

Queued independent reconstruction of all four conditional endpoint costs for
every edge, every tree and both coding scenarios. The verifier uses integer
network maximum flow rather than the producer's inside/outside recurrence.
It also reconstructs every split identity, tree score, status, per-tree count
and split-frequency denominator directly from source trees and classifications.
No producer edge-mapping helper is imported. The final proof requires all 2,000
trees; per-tree checkpoints alone do not establish completion.

A complete real bootstrap shard passed 8,392 constrained-cost checks in 6.92
seconds. Deliberately altered cost and duplicate-edge records were rejected
even with updated artifact hashes. The full run has a pre-launch 1–6-hour
allowance, four CPU workers, 8 GiB RAM, no swap and 0.1 GiB proof output. The
rough timing extrapolation is about one hour before workload variability and
summary checks. No GPU or paid service is used.

The verifier checks the producer PID plus creation time and command while
waiting, then requires its final source-bound receipt. It is recorded in
`metadata/ecology_bootstrap_edges_readback_{plan,launch}_20260927.json`;
fixture evidence is `metadata/ecology_bootstrap_flow_fixture_20260927.json`.
Full readback output is
`results/ecology/bootstrap-optimal-edge-readback-20260927-v1`.

```bash
python scripts/check_ecology_bootstrap_flow_readback.py
python scripts/readback_ecology_bootstrap_edges.py --plan metadata/ecology_bootstrap_edges_readback_plan_20260927.json
```

Production and the full independent audit remain pending. A verification rerun
recomputes the network flows; partial proof files are not silently accepted as
completed verification.

## Bootstrap production complete; full audit running

The producer exited successfully after all 2,000 trees and 4,000 coding/tree
combinations. It emitted 4,196,000 edge records and 5,006 source/coding/split
summary rows. All 2,000 compressed tree-shard hashes, their checkpoint receipts
and both summary hashes were checked against the final source-bound receipt.
Completion evidence is archived in
`metadata/ecology_bootstrap_edges_production_completed_20260927.json`.

The independent verifier automatically advanced from its dependency wait and
was observed computing with four live workers (PIDs 2167140–2167143 under the
recorded controller). Its first ten completed trees cover 83,920 constrained
cost checks. The full audit is still pending; these partial checks do not
validate the entire ensemble, and no final bootstrap transition interpretation
is yet claimed. No restart or change to pinned live code was needed.
