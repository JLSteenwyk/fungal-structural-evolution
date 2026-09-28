# Matched-model uncertainty calibration

The current conditional coefficient covariance treats the estimated variance
ratios as fixed. Neither successful optimization nor likelihood replay establishes
interval coverage when those ratios are estimated, especially at zero-component
boundaries. The completed original grid and ongoing refined grid must retain
that limitation. Polynomial ML comparisons also do not inherit calibrated
significance simply because their optimizers pass numerical checks.

## Implemented simulation primitive

`scripts/simulate_matched_working_model.py` generates responses at a fixed design:

    y = X beta + sqrt(scale) * (e + sqrt(r_bg) Z_bg u_bg
                                + sqrt(r_family) Z_family u_family
                                + sqrt(r_species) F u_species)

Each latent vector consists of independent standard normal draws. Thus the
response covariance is

    scale * (I + r_bg Z_bg Z_bg' + r_family Z_family Z_family'
               + r_species F F').

The routine uses group indexing and the existing species factor rather than
forming a dense observation covariance. It requires background groups nested
within family components, positive residual scale, finite compatible inputs
and nonnegative component ratios. It accepts an explicit NumPy random generator;
future production must record seeds and stream assignments. Return dimensions
are observations by draws. Callers should batch draws to bound memory.

The checker exercises all eight zero/nonzero component combinations with both
rank-zero and rank-two species factors. Deterministic positive/negative latent
basis fixtures reproduce the mean and covariance exactly up to arithmetic
roundoff: maximum absolute covariance error 3.55e-15 across 16 cases. Seven
invalid-input cases are rejected, and identical seeds reproduce identical
responses. These fixtures test implementation; they are not biological
replicates, a reduced sampling experiment or evidence of model adequacy.

Reproduce the checks with:

```bash
OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python scripts/check_matched_working_model_simulator.py
```

The checksum-bound result is
`metadata/matched_working_model_simulator_checks_20260928.json`.

## Remaining calibration work

Production must reconstruct each exact input design, ordered observation
identity, grouping and tree factor used in the full fit inventory. Preserve
the mapping from 144,040 unique tree/input fits to 414,720 settings; repeated
settings cannot count as independent calibration evidence. Integrate the audited
refinements first and retain all unresolved numerical cases explicitly.

Simulations used to evaluate fitted-parameter uncertainty must refit variance
components, preserve failed replicates and evaluate coverage against the known
generating coefficients. Fixed-ratio simulation alone can check the conditional
GLS calculation but cannot answer that question. Null simulations for tests
must enforce their stated null and refit both models using the same observation
set; REML objectives from different fixed-effect designs are not a valid model
comparison. Monte Carlo uncertainty and multiplicity need explicit treatment.

Before launching production, calculate resource requirements from the completed
fit timings, selected replicate counts and full unique-input inventory. No
large simulation batch has been launched or declared calibrated by this work.

Gaussian simulation under the fitted working model tests operating properties
under that model. It does not validate Gaussian residuals, constant variance,
the selected covariance form, phylogenetic assumptions or the conditioning on
observed structural data. Residual/model-adequacy diagnostics, missingness,
prediction uncertainty and biological controls remain separate requirements.

## Full-grid resource inventory

`scripts/estimate_matched_calibration_resources.py` verified all 144,040 saved
fit hashes against their manifest, checked input/tree identities and the
completed output-audit binding, and extracted per-fit elapsed times and serialized
sizes. Their summed elapsed time is 461,195 seconds (128.11 hours), with a
3.29-second median and 12.18-second maximum. These are individual elapsed
measurements under the original scheduling, not measured CPU consumption.

The scenarios below preserve all unique fits and assume unchanged fit cost,
ideal 32-worker parallelism and one generating-model condition per fit:

| Simulations per fit | Total refits | Ideal parallel days | Same-format outputs, GiB | MC standard error for 95% coverage |
|---:|---:|---:|---:|---:|
| 199 | 28,663,960 | 33.2 | 336.5 | 1.55 percentage points |
| 999 | 143,895,960 | 166.6 | 1,689.4 | 0.69 percentage points |
| 1,999 | 287,935,960 | 333.5 | 3,380.4 | 0.49 percentage points |

These are arithmetic planning scenarios, not ETAs or a selected production
design. They exclude simulation, repeated data reconstruction, audits, failed
replicate retries, extra null/alternative conditions and changed optimizer cost.
Storage scales the existing serialized fit size, not simulated alignments or
all possible auxiliary outputs. The complete scenarios include 8- and 16-worker
alternatives in `metadata/matched_calibration_resources_20260928.json`.

The original likelihood already caches invariant cross-products within a fit.
Further savings require demonstrated reuse across responses or faster equivalent
refitting, not an assumption that the original code was uncached. Simulation
responses have different fitted variance ratios, so a shared fixed covariance
alone cannot replace per-response optimization for fitted-parameter calibration.
Evaluate any batched implementation against direct likelihood and refit results
before using its timing for production estimates. The full grid remains in scope;
no smaller subset is being substituted or declared calibrated.

## Shared-covariance evaluation

`scripts/batched_matched_reml.py` evaluates multiple response columns at one
fixed variance-ratio vector, reusing the covariance solver and fixed-effect
information factorization. Each response retains its own fitted coefficients,
residual quadratic, profiled scale, conditional coefficient covariance and
REML objective. Residuals are evaluated directly rather than subtracting two
large cross-products. Invalid residual quadratics receive an explicit invalid
mask and unavailable objective/scale, without dropping other responses.

The checker compared 112 valid response evaluations against both the original
scalar evaluator and independently assembled dense covariance calculations,
covering every zero/nonzero component combination with rank-zero and rank-four
species factors. Maximum absolute objective disagreement was 7.11e-15. A zero
response stayed explicitly invalid in each batch, and four malformed input
cases were rejected.

In three single-thread implementation timings with 600 observations, 64
responses, five design columns and species rank 20, batch evaluation took
0.0038–0.0042 seconds versus 0.0273–0.0274 seconds for separate evaluations,
a 6.56–7.25-fold ratio. Both timings reuse an already constructed covariance.
They do not benchmark the existing cached optimizer, factor construction,
large-rank inputs or per-response variance optimization. Consequently they
do not revise the full-grid refitting estimate. The implementation provides
reuse at common starts/grid points and conditional-calculation checks; an
independent optimizer for each response is still necessary for refit calibration.

Reproduce with `OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python
scripts/check_batched_matched_reml.py`. Checksums and measurements are in
`metadata/batched_matched_reml_checks_20260928.json`. Existing production
scripts and running jobs were not modified.

## Independent simulation refits implemented

`scripts/refit_matched_simulation.py` connects one simulated response to the
existing analytic-gradient optimizer. Each response gets its own variance-ratio
fit, all 24 candidate records and direct candidate-likelihood checks. The
generating ratios provide an additional start and retained reference; they are
not fixed during optimization. Generating ratios must lie within the stated
optimization bounds. Random streams depend on the master seed, source-fit
identifier and replicate index, so scheduling order does not define the stream.
Each disposition records its seed entropy and response checksum.

Numerically qualified refits retain coefficient errors, fitted standard errors,
studentized errors and whether a nominal 95% t interval covers the known
generating coefficient. Those intervals are a target of calibration, not an
assumption of valid coverage. Refits needing numerical review and exceptions
remain explicit outcomes. The accounting reports lower/upper coverage fractions
obtained by treating every unresolved fit as uncovered/covered; these bounds
are not Monte Carlo confidence intervals. They prevent a success-only
denominator from hiding failures.

The end-to-end checker passed six synthetic refits spanning all-zero and
positive generating variance components, with all 144 candidate evaluations
retained and all six refits passing numerical checks. An injected optimizer
failure remained in the attempted denominator; duplicate replicate accounting
was rejected and seed/fit identifiers produced distinct streams. These are
implementation fixtures, too few to estimate coverage, and are not a biological
pilot or replacement for the full data analysis.

Reproduce with `OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python
scripts/check_matched_simulation_refitting.py`. Outputs and hashes are recorded
in `metadata/matched_simulation_refit_checks_20260928.json`. Production input
reconstruction, restartable full-grid scheduling, replicate-count selection,
Monte Carlo intervals and runtime reduction still remain. No full-grid
simulation has been launched by these checks.

## Complete simulation-input cache running

The input reconstruction stage is now running under
`fungal-matched-simulation-input-cache-20260928.service`. Its plan and verified
launch identity are in `metadata/matched_simulation_input_cache_plan_20260928.json`
and `metadata/matched_simulation_input_cache_launch_20260928.json`.

The source inventory contains 28,808 distinct inputs and 130,910,712 record
instances. Each cache entry preserves the five numeric columns, ordered record
identity hashes, nested background/family indices, species-pattern row indices,
active covariates and scale factors. Numeric and ordered identity bytes must
match the original recipe hashes. All five original phylogenetic factors are
copied once with their source hash bindings. Each serialized array is compared
exactly to its reconstructed input, including when reusing an existing cache
entry. Completed input files survive interruption without overwriting.

The resource allowance is one CPU, 16 GiB RAM, no swap and 32 GiB storage with
a 100 GiB free-disk reserve. Raw per-record arrays total about 15.61 GiB before
compression; the broader storage allowance includes files and factors. The
0.1–12-hour duration is a planning range, not a measured ETA. This stage does
not fit any models or launch simulated refits.

Complete-cache accounting and an independent likelihood replay from cached
inputs remain required before production simulation. The running cache job
is not itself evidence that those checks have passed.

## Full cache likelihood replay queued

`scripts/audit_matched_simulation_cache.py` now waits on the pinned cache
producer identity and requires inactive/success/exit0 plus complete source
receipts. It will replay all 144,040 original fixed-parameter fits from the
28,808 serialized input arrays. Checks include active covariates/scales,
phylogenetic factor indices, likelihood, coefficients, residual scale and
quadratic, conditional covariance, and conversion back to raw covariate units.
Every cached numeric and ordered-identity hash must match its source recipe.
The direct likelihood evaluator is shared with the earlier fit validation;
this is a separate input-route replay, not an independent statistical model.

Before launch, 12 available input caches across all five trees reproduced 60
original fits, with maximum objective disagreement 4.55e-12. Changing a response
value caused rejection. The exact source-bound fixture set can be reproduced
with `scripts/recheck_matched_cache_replay_fixtures.py`; evidence is recorded in
`metadata/matched_simulation_cache_replay_checks_20260928.json`.

The queued full replay uses one CPU, 8 GiB memory, no swap and 1 GiB output
allowance, with 1–24 active hours budgeted after the dependency wait. Plan and
launch identity use the `metadata/matched_simulation_cache_replay_` prefix.
The full replay has not yet passed, and no calibrated uncertainty claim follows
from this limited helper check.
