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
was unchanged and verified live by PID, creation time and command at that
snapshot. Its later terminal state is recorded below.

## Full batch terminal readback, September 28

The producer has now attempted all 156 inputs and exited with status 1:
153 nonempty inputs failed its output validator, and three inputs contained
no coded characters. All 153 inference subprocesses themselves exited zero.
The exported Newick represents the root label as comment `[N1]`; the original
validator required a node name and rejected these trees before checking the
probability rows. No producer receipts or original outputs were changed.

`scripts/audit_completed_fastml_indel_outputs.py` explicitly reads that root
comment as the root identity and checks all input/output hashes, node names,
complete position-by-node grids, probability bounds, and probabilities at
known observed tip states. This passed across all 153 nonempty inputs and
8,058,340 probability rows. Unknown tip characters are not required to have
probability zero or one. Evidence is in
`metadata/fastml_indel_terminal_readback_20260928.json`.

This resolves the output-identity parsing failure, not the likelihood
discrepancy above. Independent replay of the entire batch, resolution of
ascertainment-cache behavior, and optimization validation remain pending.
These outputs are not a qualified ancestral ensemble.

```bash
python scripts/audit_completed_fastml_indel_outputs.py \
  --output metadata/fastml_indel_terminal_readback_20260928.json
```

Use a fresh output path for a repeated audit. The command verifies the
preserved failed producer state as well as the output artifacts.
