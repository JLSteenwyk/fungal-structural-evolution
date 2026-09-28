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
discrepancy above. At this readback checkpoint, independent replay of the
entire batch remained pending; its subsequent results are below. Resolution
of ascertainment-cache behavior and optimization validation remain pending.
These outputs are not a qualified ancestral ensemble.

```bash
python scripts/audit_completed_fastml_indel_outputs.py \
  --output metadata/fastml_indel_terminal_readback_20260928.json
```

Use a fresh output path for a repeated audit. The command verifies the
preserved failed producer state as well as the output artifacts.

## Completed independent replay, September 28

The separate replay service completed successfully (inactive, exit zero).
All 153 nonempty inputs and 8,058,340 probabilities were recomputed from
serialized parameters and branch lengths; the three empty inputs remain
explicit. The analytic enumeration check also passed. Source, plan, and
per-case output hashes were verified after completion.

The largest absolute node-probability difference is 0.0001005571, in
`OG0000972-alignment-famsa-terminal_gap`. Six likelihood differences exceed
the diagnostic trigger of 0.01. Three are the OG0000294 cases tabulated
above; their discrepancies remain about 0.47, and substituting the original
tree's exclusion denominator reproduces the printed likelihood within
2.19e-6. The other three are whole-protein OG0000972 cases, with differences
0.01055–0.04759 and reported likelihoods rounded to one decimal place.
The stale-denominator diagnostic does not explain those three cases.

Across all 153 cases, 21 differences exceed half the final printed digit of
the reported likelihood. That count considers only likelihood display
precision. It does not account for rounded fitted parameters or tree edges,
so it is neither a count of confirmed software defects nor a validation
threshold. Full-precision exports or independent corrected fits are needed
to resolve these numerical differences and qualify parameter estimates.

The [complete discrepancy table](tables/fastml_complete_replay_discrepancies_20260928.tsv)
retains every case, its printed likelihood precision, probability discrepancy,
and available initial-tree-denominator comparison. Completion evidence is
`metadata/fastml_complete_replay_completed_20260928.json`.
The scripts `replay_completed_fastml_indels.py` and
`summarize_complete_fastml_replay.py` reproduce the replay and summary using
the pinned plan `metadata/fastml_complete_replay_plan_20260928.json` and new
output directories. No source predictions or failed producer receipts were
replaced. These results do not establish a qualified ancestral ensemble.

## Paired source regression, September 28

Two isolated binaries were built from the frozen source: one exports model
parameters, probabilities and tree lengths with 17-digit precision; the other
additionally refreshes the exclusion-probability cache and optimizer baseline
immediately after initial maximum-parsimony branch scaling. The original
source and binary were preserved. Patches and build recipes are tracked under
`patches/fastml-precision-cache-20260928/` and
`scripts/prepare_fastml_precision_cache_builds.py`.

Six regression runs cover the three diagnosed OG0000294 inputs under both
variants. The precision-only version retains likelihood discrepancies of
0.469812–0.473948. With the cache refresh, the largest absolute discrepancy
is 1.232e-7, below the predeclared 1e-6 regression threshold. The largest
probability discrepancy is 1.592e-8 across both variants. Increasing print
precision alone therefore does not resolve the diagnosed likelihood error;
the targeted cache refresh resolves it to the specified tolerance in these
three cases.

These checks include independent pruning/inside-outside replay and the small
enumeration check. They do not certify other cache or optimizer routes.
Full-grid paired inference and replay remain required before accepting these
variants for scientific use. The six runs are preserved in
`results/ancestral/fastml-cache-regression-20260928-v1`, with hash-checked
evidence in `metadata/fastml_cache_regression_completed_20260928.json`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/qualify_fastml_cache_regression.py \
  --output results/ancestral/fastml-cache-regression-20260928-v1
```
