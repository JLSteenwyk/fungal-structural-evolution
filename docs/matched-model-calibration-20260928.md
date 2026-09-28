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

## Monte Carlo interval accounting implemented

`scripts/matched_calibration_intervals.py` adds fixed-sample exact binomial
(Clopper–Pearson) intervals to replicate accounting. If k refits demonstrably
cover the generating coefficient and u remain unresolved among n attempts,
the lower endpoint uses k successes and the upper endpoint uses k+u successes.
This contains the exact interval for every possible resolution of the unknown
outcomes. Failed refits are not removed from n. An entirely unresolved batch
returns [0,1]. These bounds describe uncertainty about a putative resolved
procedure; they do not turn a procedure with numerical failures into a qualified
estimator. The operational rule for failures must also be reported.

`summarize_replicates` connects this calculation to the refit dispositions,
retaining duplicate-replicate rejection and separate coefficient accounting.
Fixed-sample independent simulation under one generating model is required.
The intervals are marginal Monte Carlo intervals, not coefficient confidence
intervals, not simultaneous grid-wide statements, and not valid under arbitrary
optional stopping. Production must fix the sample count or implement a separately
justified sequential method and address any intended multiplicity claims.

The checker enumerated binomial coverage at 412 sample-size/probability settings;
the minimum enumerated coverage for nominal95% intervals was95.43%. It also
checked unresolved-outcome interval containment, boundary counts, very small
alpha, six invalid-input cases and the refit-summary connection. This finite
enumeration is a numerical check, not a proof over every probability or project
coverage evidence. Reproduce with
`scripts/check_matched_calibration_intervals.py`; evidence is in
`metadata/matched_calibration_interval_checks_20260928.json`.

## Complete-refit timing running

The full simulation-to-refit path is being timed on pinned existing designs,
using minimum, median-ranked and maximum record counts among the 9,748 inputs
already cached at selection time: 148, 4,909 and 10,960 records. Each is evaluated
under all five trees, for 15 timing cases. Selection is based on availability
and size, so this is not a representative runtime sample or biological pilot.
One simulated response per case cannot estimate coverage.

Every case first reproduces the original fixed-parameter likelihood from the
cached input. It then generates a new response from the original fitted model
and independently optimizes its variance components, retaining all24 candidates
and source hashes. Original loading/replay time and simulation/refit time are
recorded separately. A first completed case took25.12 seconds and passed
numerical checks; this already demonstrates why a fast fixed-covariance
evaluation is insufficient evidence of a fast full refit.

`fungal-matched-simulation-refit-timing-20260928.service` uses one CPU,8GiB
memory,no swap,0.1GiB storage and a two-hour service runtime limit. The active
planning range is1–120minutes. Plan and launch metadata have the prefix
`matched_simulation_refit_timing_20260928`; the script is
`scripts/benchmark_matched_simulation_refits.py`. Review all timing dispositions
and failures after terminal completion before revising any full-grid estimate.

## Exact zero-species-gradient shortcut checked

Inspection of the full analytic gradient identified redundant factorization
of `I + p K` when the species variance ratio p is exactly zero. The separate
`scripts/matched_reml_gradient_fast.py` uses the identity inverse in that case,
retaining the derivative with respect to p through the cross-product and trace
terms. It also omits products multiplied by exactly zero. Positive p follows
the original calculation; there is no near-zero threshold or boundary change.
The original module remains unchanged for every active job.

`scripts/check_matched_reml_zero_shortcut.py` compares the new calculation
against the original and an independent dense residual-projector score over
81 zero/moderate/large-variance cases, including empty and duplicated species
factors. All243 derivatives passed; maximum dense-score disagreement was
7.57e-10. With600 observations and species rank250,20 fixed-point evaluations
took0.0210 seconds with the shortcut versus0.1295 seconds originally when p=0
(6.16-fold). At p=0.3 the timings were essentially unchanged.

This identifies an exact computational saving, not a6-fold full-refit speedup.
Full optimizer equivalence, boundary diagnostics and end-to-end timings still
need checking before adoption. The result is recorded in
`metadata/matched_reml_zero_shortcut_checks_20260928.json`.

## Original timing complete; paired optimizer comparison running

The original full-refit timing service reached inactive/success/exit0. All15
source-bound dispositions,24-candidate records, statuses and timing fields were
read back; every simulated fit passed its numerical checks. Timings by input
size were1.41–1.83seconds (148records),12.04–16.12seconds (4,909records), and
25.08–28.36seconds (10,960records). Completion evidence is in
`metadata/matched_simulation_refit_timing_completed_20260928.json`.

The separate `scripts/refine_matched_reml_analytic_fast.py` preserves the same
optimizer,24 candidates,bounds,tolerances and direct candidate checks, changing
only the imported gradient implementation. The queued paired comparison now
recreates each of the15 exact response hashes, reruns this optimizer and records
likelihood, coefficient, conditional covariance, status and diagnostic-flag
agreement alongside runtime. Discrepancies remain explicit; they do not silently
qualify the shortcut. No live production optimizer was replaced.

The comparison service is `fungal-matched-refit-shortcut-comparison-20260928`,
with one CPU,8GiB memory,no swap,0.1GiB output and1–120active minutes budgeted.
Plan and launch metadata use `matched_refit_shortcut_comparison_`. Results are
pending and these size-selected cases cannot establish global equivalence or
a representative whole-project speedup.

## Paired optimizer comparison completed

The comparison reached inactive/success/exit0, and all15 paired output hashes,
candidate counts, estimate comparisons, diagnostic flags and timing arithmetic
were read back. All checked agreements passed. The largest absolute fitted
objective difference was1.82e-12. Original total measured time was210.94seconds
versus196.03seconds with the shortcut, a7.07% reduction over these cases.
Individual speed ratios ranged from0.99 to2.52. There is only one timing per
response, so host-load variation remains and no statistical runtime claim is
made. Evidence is in
`metadata/matched_refit_shortcut_comparison_completed_20260928.json`.

The shortcut is a modest implementation improvement on this selected set. It
does not justify applying the earlier6-fold fixed-point improvement to the
full calibration estimate, and it does not resolve the very large refit count.
Existing active production fits were not replaced or restarted.

## Analytic uncertainty feasibility and full boundary census

The downloaded [Halekoh & Højsgaard (2014) article](https://www.jstatsoft.org/article/view/v059i09),
Appendix A, expresses marginal covariance as a linear combination of known
matrices and supplies information, fixed-effect covariance adjustment and
F-moment calculations. This suggests a possible analytic alternative to many
simulation refits, subject to independent validation. It does not establish
valid coverage for this project.

Our model can be written in absolute variance coordinates as
V = v0 I + vb Zb Zb' + vf Zf Zf' + vs F F', where v0 is the profiled scale
and each other variance is that scale times its fitted ratio. Thus there are
four known covariance kernels. Any information calculation must include the
residual variance; the three-ratio profiled optimizer Hessian alone is not the
four-parameter covariance required for this approach. This is our algebraic
mapping, not a claim that the existing software directly supports our fitter.

Before choosing a method, `scripts/inventory_matched_variance_boundaries.py`
checked all144,040 original unique tree fits against the hashed export. There
are61,932 fits (43.00%) with at least one exactly zero variance component:
58,220 have zero family variance,24,312 zero species variance, and20,600 both.
Background variance is positive in every fit. All658 numerical-review fits are
retained; the ongoing refinement overlay is not part of this census. No
near-zero threshold was applied. Evidence and per-tree accounting are in
`metadata/matched_variance_boundaries_20260928.json` and the bound output files.

These counts make boundary handling a central validation requirement. Zero
estimates do not automatically invalidate a fit or prove that an approximation
fails. Conversely, the82,108 fits with positive components are not automatically
identifiable or calibrated. Information rank, conditioning, fixed-effect
estimability and coverage remain to be checked. Do not silently drop boundary
components, exclude these fits from denominators, or report an approximation
as calibrated because it returns finite numbers.

One dense float64 observation covariance at the largest input would occupy
0.8955GiB, before temporary matrices. A full-grid implementation should exploit
the nested groups and low-rank species factor rather than repeatedly construct
dense kernels. Next implementation work should validate exact information and
covariance contractions against independent dense calculations, including
rank-deficient and boundary cases, then compare approximation behavior with
refitted simulations. No full calibration replicate count or analytic method
has yet been selected; the entire original analysis grid remains in scope.

## Exact information reference implemented

`scripts/matched_covariance_information.py` now computes the four-component
expected REML information and fixed-effect covariance contractions in absolute
variance coordinates. With T=V^-1 X, Phi=(X' V^-1 X)^-1 and
P=V^-1-T Phi T', the information is Iij=tr(P Gi P Gj)/2.
For Gi=Ui Ui', its entries are accumulated as
tr((P Ui)' Gj (P Ui))/2 over blocks of columns of Ui. Identity, group-incidence
and species-factor kernels are applied without storing observation-sized dense
covariance matrices. The residual identity contribution still requires
quadratic work; this is a memory-bounded exact reference, not yet a scalable
production implementation for144,040 fits.

The same routine returns -T' Gi T and T' Gi V^-1 Gj T, together with Phi,
for future covariance-adjustment calculations. All derivative kernels remain
present when their fitted variances are zero. It reports numerical rank from
diagonally normalized information and its eigenvalues, with an explicit
floating-point tolerance. It does not invert information, remove confounded
components, calculate degrees of freedom or issue intervals.

`scripts/check_matched_covariance_information.py` passed48 independent dense
comparisons spanning all eight zero/nonzero random-component combinations,
two overall scales, ordinary/duplicated/empty species factors. Maximum absolute
disagreement across all returned contractions was3.73e-9 under the recorded
relative/absolute checks. A species kernel identical to the family kernel
correctly gives numerical information rank3. Changing block partitions gave
equivalent information, and seven invalid inputs were rejected. Evidence is
in `metadata/matched_covariance_information_checks_20260928.json`.

These are numerical checks of the formulas, not coverage evidence. Boundary
sampling behavior, conditioning on real inputs, an efficient full-grid path,
covariance adjustment, reference-distribution calculations and simulation
validation remain pending. Active production scripts remain unchanged.

## Real-input information timing running

The exact reference is now running on all15 existing size/tree timing designs
(148,4,909,10,960 records, each under five trees). These are the previously
selected performance cases, not a new biological subset or a coverage sample.
Every case replays the original fitted likelihood and coefficients from the
serialized simulation input, then evaluates information at the original
absolute variance estimates. Its conditional coefficient covariance must also
match the original fit. Full matrices, normalized information eigenvalues,
numerical rank, zero-component indicators, timings and source hashes are saved.
Serialized arrays are read back exactly before each case is recorded.

The service `fungal-matched-information-timing-20260928.service` has one CPU,
8GiB memory,no swap and a two-hour runtime limit. Planning allows1–120minutes
and0.1GiB outputs; one n-by-64 work buffer at the largest input is5.35MiB,
with additional factor and solve workspaces. Launch PID, creation time and
command were verified, and CPU accounting confirmed active work. Plan and
launch metadata use `matched_information_timing_20260928`; outputs are in
`results/model_validation/matched-information-timing-20260928-v1`.

Completion and artifact auditing remain pending. These original estimates do
not incorporate the separate658-fit refinement overlay. A full-rank result
will not by itself establish regularity or valid uncertainty at a boundary.

## Residual information trace shortcut checked

The separate `scripts/matched_covariance_information_fast.py` avoids the
observation-identity block pass. It computes tr(P squared) using the nested
base inverse and the species/fixed-effect low-rank downdates. Cross-information
with the residual component is accumulated as squared norms of P Ui during
the existing group/factor passes. All four components and derivative kernels
remain present. This changes computation, not model specification.

For a downdate A-H M H', the squared trace is evaluated as
tr(A squared)-2 tr(M H' A H)+tr((M H' H) squared). The nested base contribution
uses background counts and family sums. If the final residual trace is
nonfinite or no larger than1e-8 times the sum of absolute trace terms, the
routine falls back to the original blocked reference for the entire result.
This conservative cancellation trigger is a numerical safeguard, not a
statistical threshold or a guarantee against every conditioning problem.

`scripts/check_matched_covariance_information_fast.py` compared162 synthetic
cases against both the original blocked information and dense calculations,
covering zero,moderate and10,000-fold variance ratios. Maximum information
disagreement with the blocked reference was3.32e-9. Across all contractions,
maximum dense absolute disagreement was1.20e-5 and maximum error scaled by
1+absolute reference value was2.34e-7; the large-ratio checks use explicit
rtol/atol1e-6. An adversarial full-rank species factor triggered cancellation
fallback and exactly reproduced the reference. Confounding, block-partition
equivalence and invalid-input checks also passed.

One600-record,242-factor synthetic timing took0.102seconds for the reference
and0.065seconds for the shortcut. This single timing is not a representative
full-grid speed estimate. Evidence is in
`metadata/matched_covariance_information_fast_checks_20260928.json`. Group
passes still require substantial work; real-input equivalence and timing are
pending. The active original15-case benchmark and its pinned code are unchanged.

## All real-input information comparisons passed

The original benchmark reached inactive/success/exit0. Its15 source-bound
outputs and artifact hashes were read back before running the shortcut on the
same designs and variance estimates. `scripts/compare_matched_information_shortcut.py`
completed with exit0; `scripts/audit_matched_information_comparison.py` then
checked every saved paired array, information eigenvalue/rank and timing sum.
All15 comparisons passed at rtol1e-7/atol1e-8. Maximum absolute information
disagreement was2.91e-10. Evidence is in
`metadata/matched_information_shortcut_comparison_20260928.json` and
`metadata/matched_information_shortcut_readback_20260928.json`.

Measured total time fell from212.02 to92.69seconds (56.28% reduction). Large
cases took14.50–15.56seconds with the shortcut versus33.92–35.77seconds for
the reference; medium cases took3.34–3.75seconds. Small cases were slightly
slower with the shortcut. Each method was timed once per case under current
host load; these selected designs do not provide a full-grid runtime estimate.
This remains an information calculation at fitted parameters, not a refit or
coverage simulation.

All15 information matrices have numerical rank4. Ten cases have at least one
zero fitted variance component. Their smallest diagonally normalized
information eigenvalue is at least0.2417. Thus boundary estimates and a
full-rank information matrix coexist in these inputs; rank alone does not
qualify the boundary sampling approximation. The next statistical work remains
covariance adjustment/reference-distribution validation and simulation-based
assessment, with all original analysis settings retained.

Separately, the full input-cache service reached inactive/success/exit0 and
produced its28,808-input completion receipt. The independent five-tree
likelihood replay is active (1,200 inputs observed at13:38 EDT); its full
qualification remains pending. No production input/model scope was reduced.

## Candidate covariance adjustment checked against author code

`scripts/matched_kr_covariance.py` implements the linear-kernel covariance
adjustment from the stored information and contractions. The information
matrix here is Iij=tr(P Gi P Gj)/2, so the variance-parameter covariance is
I^-1, not twice that inverse. This convention was checked against the
[authors' implementation](https://github.com/hojsgaard/pbkrtest/blob/bb7c4f9a91068b090aa9cce3c3e1d6e2de67399e/R/KR_vcovAdj.R),
which inverts twice-information and applies the corresponding factor of two.

The local comparison executes the unmodified author covariance routine and
index helper at pinned commit `bb7c4f9a91068b090aa9cce3c3e1d6e2de67399e`,
using R4.3.3,Matrix1.6.5 and MASS7.3.60.0.1. It does not install pbkrtest or
claim an independent parameter fit. Source files, retrieval URLs and hashes
are bound by `metadata/pbkrtest_covariance_reference_source_20260928.json`.

All16 synthetic comparisons (eight zero/nonzero component combinations at two
scales) passed for both adjusted fixed-effect covariance and variance-parameter
covariance, with maximum absolute difference7.99e-15. The comparison uses the
same designs and absolute variances in both implementations, with independently
constructed dense covariance kernels for the R routine. Reproduce using
`scripts/check_matched_kr_covariance.py`; full input/output hashes and R versions
are retained under `results/model_validation/matched-kr-covariance-checks-20260928-v1`.

The project implementation returns an explicit review status when information
has nonpositive diagonal entries or its normalized minimum eigenvalue is no
greater than1e-10 times its maximum. It does not use a generalized inverse;
this differs deliberately from the upstream fallback. A confounded-kernel
test verifies this behavior. Adjusted covariance must also be finite and
positive definite. These are numerical eligibility checks only. Every returned
adjustment remains a candidate pending statistical validation, and exact-zero
variance flags are retained. Degrees of freedom, reference-distribution
validation, simulation coverage, multiplicity and model adequacy are unresolved.

## Scalar moment calculation checked

`scripts/matched_kr_scalar.py` adds rank-one contrasts, including the matched
intercept, to the candidate covariance adjustment. In the rank-one moment
formulas A1=A2, so denominator degrees of freedom simplify to2/A2 and the
F scaling factor is1. The implementation uses this algebraic simplification
and the adjusted contrast variance to compute candidate t limits. Nonpositive
or nonfinite moments and inherited covariance-review dispositions remain
explicit. Joint, higher-rank tests are not implemented by this scalar routine.

`scripts/check_matched_kr_scalar.py` compares two contrasts on each of the16
existing synthetic covariance designs with the unmodified author `.KR_adjust`
function at the previously pinned commit. All32 comparisons of denominator
degrees of freedom, F scaling, adjusted variance and candidate interval limits
passed; maximum absolute disagreement was5.26e-13. Four invalid contrasts or
levels were rejected. Source provenance is in
`metadata/pbkrtest_scalar_reference_source_20260928.json`, with check results
in `metadata/matched_kr_scalar_checks_20260928.json`.

The retained v1 fixture run checked the same numerical formulas but incorrectly
filled its boundary-flag metadata with false values. The v2 fixture corrects
the flags from the original generating variance grid and verifies they survive
into candidate results; no model parameters or numeric intervals changed.
The v2 output is `results/model_validation/matched-kr-scalar-checks-20260928-v2`.
Formula agreement does not demonstrate nominal coverage. Actual refitted
simulation, boundary behavior, full-grid qualification and multiplicity remain
required before biological inference.

## Refitted-simulation interval path connected

`scripts/calibrate_matched_kr.py` evaluates candidate intervals at each
simulation's newly fitted variance components and scale. It first checks the
recomputed conditional covariance against the refit. The generating coefficients
are used only to assess whether each candidate interval contains its known
target. Failed or numerically unresolved refits are retained without intervals;
information/interval failures remain unresolved by coefficient. Evaluation
exceptions retain the attempt and record the error rather than reducing the
denominator.

The summary requires distinct replicates of one fit identifier with the same
coefficient truth and interval level. It reports the existing fixed-size exact
binomial interval envelope with unresolved outcomes retained. The calling
simulation plan must additionally bind the generating covariance/design and
prespecify the replicate count; identifiers alone cannot prove a shared
generating model. These summaries are marginal, not multiplicity adjustments,
and are not valid under optional stopping.

`scripts/check_matched_kr_refit_accounting.py` replayed all six saved synthetic
responses by seed and SHA256, then independently reproduced their fitted
likelihood, coefficients, scale and covariance. All12 candidate intervals were
computed at refitted variances. Injected calculation errors, both upstream
review statuses, coefficient-specific failures and three invalid pooling cases
were checked. All attempts remained accounted for. Evidence is in
`metadata/matched_kr_refit_accounting_checks_20260928.json`.

This completes the implementation connection, not coverage calibration. The
six fixtures are insufficient for a coverage claim; prespecified simulation
experiments, full-grid numerical qualification and model adequacy still remain.

## Prespecified synthetic coverage experiment running

`scripts/prepare_matched_kr_simulations.py` freezes two Gaussian designs
(48records,16backgrounds,4families,4species-factor columns,2fixed coefficients;
240records,60backgrounds,12families,16factor columns,5coefficients). Each has
all eight combinations of zero or0.7 background/family/species variance ratios
at residual scale1.3. There are999 seeded independent responses per configuration,
for15,984 full24-candidate refits. This is method validation, not a biological
pilot or a replacement for the full fungal sampling and analysis grid.

Replicate count and designs were fixed before launch. At coverage0.95,999
replicates give binomial Monte Carlo standard error about0.00690. All attempted
refits and interval dispositions remain in the fixed denominator. End-of-run
summaries use marginal exact binomial interval envelopes, including unresolved
outcomes. Cross-configuration comparisons are descriptive; no optional stopping,
joint significance claim or universal real-input coverage follows from this run.

`scripts/run_matched_kr_simulations.py` verifies source/design hashes, locks the
output, binds every replicate to the plan and stores the entire refit plus
interval result. Resume requires matching plan, payload checksum and replicate
identity; incomplete payload/checksum pairs require review. The controller
checks disk reserve and observed output size every100 dispositions. A completion
receipt requires all planned draws and summaries, with a separate final audit
still necessary.

The service `fungal-matched-kr-simulations-20260928.service` was launched with
four CPU workers,8GiB memory,no swap and a24-hour service ceiling; no GPUs or
paid resources. Prelaunch planning allowed1–24hours and4GiB output with100GiB
free reserve. The estimate deliberately spans0.5–20seconds per refit; initial
small-design processing is faster, but does not establish later-case throughput.
Controller and four worker identities were verified live, and598 saved
dispositions were observed at the first inspection. No interim coverage-based
decision was made. Plans and launch records have the prefix
`metadata/matched_kr_simulations_`; results are under
`results/model_validation/matched-kr-simulations-20260928-v1`.
