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

The [independent spectral implementation](../scripts/independent_shared_entity_likelihood.py)
replays candidates using component eigen decompositions, species SVD whitening
and a separate fixed-design SVD. It imports neither production likelihood nor
covariance code. Inverse residuals are checked against the supplied covariance;
unresolved component eigenvalues are review errors without clipping or jitter.
Entity traces stream at most 32 columns. REML traces use squared norms after
orthogonal fixed-effect removal. Coefficient/covariance comparisons use balanced
design units to reject corrupted tiny-unit exports. Objective, scale, variance
components, parameter identity, kernel normalization, boundary scores and
numerical-only inference flags are checked separately.

Local curvature screens evaluate independent scores at two finite-difference
resolutions with central or bounded one-sided stencils. Strictly positive lower
boundary scores permit removing those directions from the critical cone;
weak lower boundaries remain in the tested space. Positive curvature is
required on the whole remaining space, a sufficient, conservative constrained
minimum condition. Flat, negative, asymmetric and unresolved curvature remains
review. The stability envelope is empirical, not a rigorous floating-point
bound. This does not prove a global optimum or calibrate variance uncertainty.

The [final independent validation](../metadata/independent_shared_entity_likelihood_validation_20261002_v5.json)
passed 48 dense ML/REML cases, permutations, exact zero variance/rank, streamed
two-column batches and extreme design units. Maximum objective/gradient
discrepancies were 2.14e-14/3.56e-14; maximum inverse relative residual was
6.41e-15. Signed 80-digit objective/gradient discrepancies were at most
4.12e-11/2.28e-17. Both optimizer replays and their local curvature checks passed;
dense curvature discrepancy was at most 7.29e-10. All 15 altered candidate
exports, seven invalid inputs and six analytic critical-cone stress fixtures
behaved as required. These fixtures are not production fits, a biological pilot
or accepted effects. Earlier unlaunched validation versions remain preserved.

The [prospective full fitting inventory](../metadata/full_shared_entity_fit_resource_inventory_20261002_v3.json)
preserves all 622,080 original settings, both outcomes/modes, five trees and
ML/REML. Original fit-input sharing permits at most 5,208,000 unique candidates
and 12,441,600 setting links. Covariance review rows remain linked even when
no optimizer runs. Three starts with 500 search evaluations and each final
replay give a ceiling of 7,827,624,000 evaluations, not expected iterations or
time. Six-ratio candidate replay/curvature uses at most 26 independent
evaluations; an independent global-optimization budget remains unspecified.
Conservative uncompressed record allowances total 488.7 GiB, with proposed
1 TiB output/scratch reserve and two CPU/32 GiB/no swap/BLAS1 per worker.
These are prospective budgets, not launched limits or measured peak memory.
Eligible counts, full-scope timing and restartable export/readback contracts
still require closed design and covariance qualification. No production fitting
job was launched by this inventory.

The [independent optimizer](../scripts/independent_shared_entity_optimizer.py)
uses spectral scores with SLSQP, distinct from production Cholesky/L-BFGS-B.
Three deterministic starts must independently meet boundary KKT conditions,
avoid the upper cap and agree with the supplied candidate at an explicit
absolute likelihood tolerance. All failures and disagreements remain review.
[Six software comparisons](../metadata/independent_shared_entity_optimizer_validation_20261002.json)
checked every finite independent objective against a dense inverse: maximum
discrepancy 2.85e-14. Five cases met every independent-start requirement; one
retained numerical review. Budget exhaustion and an altered candidate were
rejected. Bounded multistart agreement does not prove a global optimum.

The [full fitting producer](../scripts/prepare_full_shared_entity_fits.py) and
[independent reader](../scripts/readback_full_shared_entity_fits.py) now connect
these algorithms to the entire original source grid. Candidate identities bind
closed qualification, original design/response/audit identities, methods,
optimizer/audit settings and all backend pins. Response/design hashes and
original ordered cohorts are rechecked before fitting. Only qualified bases
are optimized. All excluded, constant, dependent and failed rows remain linked.
Each cohort's compressed candidates and receipt are finalized atomically;
restart verifies and reuses completed cohorts without rewriting them. SQLite
links every original setting to its exact shared candidate. Completed producers
and readers refuse reruns, including a different reader output name.

The reader reconstructs every original candidate and setting link, independently
checks coefficients/scale/components and every finite start, validates seeds,
selected parameters and spread tolerances, and deterministically reproduces
claimed production failures. Converged candidates additionally receive local
curvature and independent SLSQP checks. Numerical reviews cannot be promoted.
Reader scratch can be rebuilt only under its exclusive lock and unchanged
source state, with a complete independent replay. Numerical output acceptance
still requires final full hashes and both actual original process journals.

[Complete software contracts](../metadata/full_shared_entity_fit_contract_validation_20261002_v3.json)
passed the 24-case/six-cohort/720-setting grid: 7,200 candidate rows and 14,400
setting links, all 12 rehashed altered-grid exports and seven altered numerical
exports, checkpoint reuse and completed-restart refusal. Separate numerical
fixtures exercise independent candidate/curvature/search qualification,
budget-failure reproduction and the full five-entity-plus-species recipe's
strict optimization review. The early SLSQP checker expectation was corrected
to retain either iteration-limit review or explicit evaluation-budget failure;
neither outcome permits numerical acceptance. Earlier unlaunched full-grid
receipts remain preserved. These are software fixtures, not a biological pilot.

The [71-pin draft production plan](../metadata/full_shared_entity_fit_draft_plan_20261002.json)
retains the full 622,080-setting grid. Its [updated resource inventory](../metadata/full_shared_entity_fit_resource_inventory_20261002_v4.json)
includes independent three-start searches, start replays and local curvature:
at most 1,534 independent-reader evaluations per candidate, including the
production-failure replay alternative. These are evaluation ceilings, not
expected work or runtime. The plan is neither launched nor queued. Closed
qualification, full-scope timing and actual resource installation are still
required; nonuniform/control variants and inferential calibration remain open.

The [full-stage environment](../environments/full-shared-entity-fits-20261002.yml)
pins the locally validated Python, NumPy, SciPy, Arrow, psutil and mpmath
versions. The smaller numerical likelihood environment alone does not supply
the full source-table/monitoring dependencies.
