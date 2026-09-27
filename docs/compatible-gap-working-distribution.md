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

## Complete candidate-set application

All936 candidate dispositions were evaluated:918 nonempty models and18
empty cases. Every nonempty case has positive compatible mass. Across77,742
character marginals,4,714 MAP states change; the maximum marginal change is
0.999921. Minimum log retained mass is -250.6935 (approximately1.3e-109).
Thus compatibility conditioning can radically alter this working
distribution. Feasibility does not qualify it as the final evolutionary
model or demonstrate accurate ancestral reconstruction.

The audited input probabilities contained65 floating-point values slightly
above1 (maximum1+1.3e-15). Only excursions outside[0,1] within1e-12 were
snapped to the boundary, explicitly counted, and retained beside raw values.
No probabilities inside the interval were regularized. Every conditional
marginal exclusion sum and every selected MAP configuration was checked.
Outputs: `results/ancestral/compatible-gap-application-20260927-v1`.
Original probabilities remain unchanged. Further work needs a biologically
justified joint indel treatment and validation, not acceptance of the
conditioned surrogate merely because it yields legal sequences.
