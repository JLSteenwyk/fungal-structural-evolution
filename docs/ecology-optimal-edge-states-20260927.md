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
