# Shared-entity likelihood and variance optimization

The new backend evaluates profile ML and REML and their analytic variance-ratio
gradients for explicit shared-entity and species kernels. It leaves the launched
covariance operator and full qualification jobs unchanged. Production variance
fitting is not launched or accepted by these software checks.

For supplied covariance `V`, fixed design `X`, GLS residual `r`,
`P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1`, `alpha=V^-1 r` and
`q=r' V^-1 r`, the REML score for kernel `K` is
`0.5 * [trace(P K) - (n-p)/q * alpha' K alpha]`.
Profile ML uses `trace(V^-1 K)` and `n/q`. Overall scale is profiled at
`q/(n-p)` for REML and `q/n` for ML. Residual quadratics are calculated directly;
perfect or numerically unresolved fits are rejected.

The backend balances active columns and uses QR coordinates for fixed effects,
then restores original-unit coefficients and their conditional covariance.
The fixed-design determinant Jacobian is retained for REML. Unrepresentable
coefficient covariance is an explicit error, including numerical underflow.
Conditional coefficient covariance does not propagate estimated variance,
prediction, phylogenetic or ancestral uncertainty. REML objectives across
different fixed designs must not be used as ordinary model-comparison
likelihoods; ML comparisons still require appropriate calibration.

Entity traces use component-local inverse blocks, with the species low-rank
correction, and explicit fixed-effect subtraction. Species traces use checked
covariance solves. Component kernels are cached for repeated evaluations.
No global dense covariance or residual projector is constructed.

The optimizer first requires independent raw and residual-space covariance
bases, evaluated against the original active design. Kernel-diagonal norms
scale variance coordinates. `ratio = expm1(u)/kernel_norm`, `u>=0`, permits
exact zero variance; logarithmic positivity floors are not imposed. Three
deterministic starts retain every evaluation/iteration/failure disposition.
Final gradient replay and an explicit boundary KKT check supplement optimizer
success flags. An upper scaled-variance cap is a computational review boundary,
not evidence of a finite optimum. Failed starts, disagreement, upper-bound
solutions and budget exhaustion remain review states. Even a converged
candidate is pending independent optimization, curvature and numerical audits.

Current [likelihood validation](../metadata/shared_entity_likelihood_validation_20261002_v3.json)
passed 36 dense ML/REML comparisons, fourth-order forward/central numerical
derivatives, row permutations, exact variance-zero boundaries, zero-rank species
factors and extreme fixed-design units. Maximum analytic-gradient discrepancy
was 7.11e-15; maximum numerical-derivative discrepancy was 1.03e-9. Nine invalid
inputs were rejected. The [80-digit signed check](../metadata/shared_entity_likelihood_precision_validation_20261002_v2.json)
used strong entity/species variances 10,000/3,000; maximum score discrepancy
was 1.90e-16 and objective discrepancy 3.30e-12.

[Optimizer contracts](../metadata/shared_entity_optimizer_validation_20261002_v2.json)
passed all six independent dense Powell comparisons, with objective
discrepancy at most 4.27e-10. Every ordinary test candidate passed all three
starts but remains pending audit. Cap-limited solutions, exhausted budgets,
unidentified target/residual terms and misleading narrow-box optimizer success
were retained as review/failure states. These are synthetic numerical software
fixtures, not a biological pilot or production inference. Earlier unlaunched
receipts are preserved; linked versions bind the current implementation.

The [environment](../environments/shared-entity-likelihood-20261002.yml) records
exact package versions. Each checker requires a new exclusive output path.
Full production timing/resource estimates, qualified cohort/mode/tree recipes,
independent full-scope fitting audits, weight/control variants and uncertainty
calibration remain required. All eight scientific aims remain incomplete;
structural GPU prediction remains paused.
