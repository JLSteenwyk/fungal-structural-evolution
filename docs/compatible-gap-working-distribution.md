# Compatible exact-gap working distribution

The independently inferred SIC marginals can violate necessary exclusion
constraints. A compatibility calculation cannot preserve incompatible
marginals exactly. We retain them as the original model outputs.

`scripts/compatible_gap_distribution.py` defines an additional, explicit
working distribution: start with product-Bernoulli weights from the supplied
gap probabilities and condition on nonoverlapping, nonadjacent exact gap
runs. Dynamic programming sums the weights of every compatible configuration,
computes the new marginals and a maximum-weight configuration, and supports
backward sampling. The retained mass and marginal changes must be reported.

The recurrence orders intervals by end position. Either omit the current
interval, or include it while omitting all conflicting preceding intervals
and continuing at the last compatible prefix. Calculations use log weights;
zero and one input probabilities are retained exactly. Zero compatible mass
is an explicit infeasible case, not repaired by silently clipping probabilities.

Exhaustive enumeration validates partition masses, marginals and MAP weights
for84 small cases, including adjacency, overlap, nesting, empty inputs,
deterministic states and impossible support. Random draws are checked for
compatibility; the current tests do not check sampling frequencies.

This distribution is a computational sensitivity model, **not a fitted joint
evolutionary indel process**. It assumes product weights before conditioning,
changes the supplied marginals, and only permits gap intervals represented
in the observed coding. It does not establish historical event counts or
support unobserved ancestral gap intervals. Its relationship to unknown-state
coding and alternative alignments remains an uncertainty source. Applying
it to the complete candidate set, measuring its effect, and evaluating a
justified joint evolutionary alternative remain necessary before treating
resulting sequence ensembles as qualified ancestral reconstructions.
