# Independent FastML replay: unresolved likelihood discrepancy

`scripts/replay_fastml_indel_probabilities.py` independently implements a
two-state continuous-time model, gamma category means, scaled pruning and
inside/outside marginal probabilities. A small tree with unknown tip states
passes independent enumeration of all node assignments for both likelihoods
and probabilities to 1e-12. The production replay uses serialized parameters
and branch lengths, which have limited precision; it does not claim a
full-precision match or convergence.

In the first 81 nonempty completed inputs, the maximum node-probability
difference is approximately 6.19e-5. Three OG0000294 one-character fits have
likelihood differences much larger than their printed precision:

| Encoding | Replay minus reported log likelihood |
| --- | ---: |
| envelope/FAMSA, terminal gap | -0.469816 |
| envelope/MAFFT, terminal gap | -0.473947 |
| envelope/MAFFT, terminal unknown | -0.473947 |

Using the final tree for the observed-pattern numerator but the original
input tree for the all-zero exclusion denominator reproduces these reported
likelihoods within approximately 2.2e-6. This is evidence consistent with a
stale ascertainment denominator after initial tree scaling. It does not by
itself prove every route through the optimizer is affected.

The inspected FastML source supports this diagnosis:
`gainLoss::multipleAllBranchesByFactorAtStartByMaxParsimonyCost` rescales the
tree without refreshing `unObservableData`. The likelihood routine consumes
the cached exclusion probability. These cases retain their initial gamma,
gain and loss parameters. Other optimizer routes may refresh the cache.

Keep the original outputs and explicitly withhold likelihood validation for
these cases. A corrected inference run must recompute the exclusion
probability whenever parameters or branch lengths change. Comparison with
the original method, missing-mask conditional sensitivity, optimization
checks and complete independent replay remain pending. Internal-node
probability agreement alone cannot validate fitted model parameters.

The named snapshot in `metadata/fastml_indel_probability_replay_snapshot_20260927.json`
contains exact discrepancies and source receipt hashes. The full producer
remains unchanged and was verified live by PID, creation time and command.
