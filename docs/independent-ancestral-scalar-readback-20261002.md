# Separate ancestral scalar numerical readback

This stage recalculates the existing scalar and candidate-length diagnostics
for the complete recovered-attempt overlay. It retains **405 original quartets**,
including **403 complete and two unresolved**, at both original burn-in cutoffs.
All 37 scalar variables and four candidate-node lengths remain in scope:
**33,046 variable/cutoff rows**. The original models, inputs, seeds, priors,
thresholds, failed attempts and published diagnostic outputs remain unchanged.
The stage completed its serialized readback and provenance closure on October 2.
This numerical check does not establish a qualified ancestral posterior.

## Numerical method and review rules

The [separate implementation](../scripts/independent_ancestral_scalar_diagnostics_v2.py)
calculates rank-normalized/folded split R-hat, bulk ESS, tail ESS and mean Monte
Carlo standard error. Rank/folded diagnostics and quantile-localized ESS follow
[Vehtari and colleagues](https://arxiv.org/abs/1903.08008). The
[Stan reference](https://mc-stan.org/docs/reference-manual/analysis.html)
describes the between/within-chain variance and autocorrelation estimators.
The implementation imports no ArviZ or project estimator. Direct lag dot
products replace the original FFT autocovariances; paired initial-positive
sequences and cumulative minima implement the monotone truncation. The declared
finite-lag refinement, antithetical ESS allowance and
[mean MCSE convention](https://python.arviz.org/en/v0.22.0/api/generated/arviz.mcse.html)
target the locked ArviZ 0.22.0 reports rather than a silently changed estimator.

Four chains are split into equal halves, dropping the central observation for
odd retained lengths. Pooled tied ranks use averaged ranks. Folded ranks use
absolute departures from the pooled median. Tail indicators use the pooled
5th/95th percentiles before splitting. The original screen still requires
R-hat strictly below 1.01 and both ESS values at least 400. A constant original
chain retains the original no-metric review disposition. Nonfinite and
insufficient traces remain explicit. The thresholds are not adjusted to make
the saved 75/50 length-state draws pass.

Exactly constant split trajectories have zero mathematical within-chain
variance. Floating summation can return an undefined or an extremely large
finite R-hat. Such a trace retains `singular_metric_requires_review`; matching
large values is not counted as successful R-hat numerical comparison. All other
defined values use absolute plus relative tolerance `1e-8`. Any unexplained
value, finite/undefined or screen-disposition mismatch fails the stage. This
distinguishes defined numerical agreement from agreement on an unresolved review.

The locked ordered maximum uses the bulk diagnostic when a globally constant
folded score contributes no defined scale diagnostic. This convention is kept
explicitly for report compatibility. None of these scalar rules verifies joint
alignment/state mixing, root uncertainty or biological adequacy.

## Complete source and serialized checks

The [source loader](../scripts/independent_baliphy_scalar_sources.py) requires
closed recovered diagnostics and their full 50,019-binding archive. It verifies
the original 1,620 chain identities, unique four-chain groups/seeds, selected
whole-attempt dispositions and scalar native sample-audit/log links. Every
report manifest is checked against its exact cutoff, variable set, expected
iteration schedule and original log hashes. Raw TSV parsing preserves every
declared variable, draw and chain; failed groups are retained without diagnostic
substitutes. Length traces inherit their prior closed native/state projection
checks; this stage does not independently reimplement native alignment parsing.

The [producer](../scripts/prepare_independent_baliphy_scalar_readback_v2.py) writes
atomic immutable quartet checkpoints containing each original and separately
calculated diagnostic, defined errors and unresolved numeric flags. On resume
it recalculates complete checkpoint contents before accepting unchanged files.
The [serialized reader](../scripts/readback_independent_baliphy_scalar_readback_v2.py)
reconstructs every row again from raw traces, rejects foreign/missing/duplicated
rows and reproduces every count partition against closed original summaries.
It uses the same separate estimator and is not a third independent numerical
implementation. Full source/artifact hashing and both actual original completion
resource journals gate final accounting closure.

## Software validation and resources

[Numerical fixtures](../metadata/independent_ancestral_scalar_validation_20261002_v3.json)
passed 107 cases covering ordinary and odd lengths, positive/negative serial
correlation, location/scale mismatch, heavy tails, categorical/binary ties,
constant and nonfinite chains. Seven exactly constant-half cases retain explicit
singular R-hat review. Maximum defined R-hat error was `2.23e-16`, bulk ESS
`6.19e-11`, tail ESS `9.10e-13` and mean MCSE `4.17e-16`. Three invalid shapes
and an altered numerical export were rejected.

The [complete workflow fixture](../metadata/independent_baliphy_scalar_grid_validation_20261002_v2.json)
passed all 405 groups and 33,046 rows, unchanged interruption-checkpoint reuse
and completed-producer/alternate-reader restart refusal. Seven rehashed false
exports were rejected, covering omitted/duplicate variables, changed metrics
or source status, foreign reports, scientific promotion and hidden groups.
Its constant traces and upstream journals are synthetic; the 107-case oracle
tests provide the separate nonconstant numerical coverage. Neither fixture is
a biological pilot or a production result.

Before launch the [plan](../metadata/independent_baliphy_scalar_plan_20261002_v2.json)
fixed two CPU equivalents, 32 GiB, no swap, one BLAS thread, a 4 GiB output
allowance and 104 GiB minimum free disk. The source archive represents 60.65 GiB
per full hash pass; multiple passes are required. Raw traces are cached per
quartet, with at most 4×1,001×37 scalar cells. The 0.5–8 hour planning range per
stage is uncalibrated and excludes scientific sampling; it is not a project or
convergence ETA. No GPU, new sampling job or paid infrastructure is used.

## Full production result

The [completion record](../metadata/independent_baliphy_scalar_completed_20261002_v2.json)
retains all 405 groups and 33,046 rows. Defined metrics were compared on 30,478
rows; the remaining 2,568 retain the original constant-chain disposition without
defined metrics. No production row required the special zero-split-variance
review. Maximum absolute discrepancies were `4.45e-16` for R-hat, `5.30e-11`
for bulk ESS, `5.96e-11` for tail ESS and `5.69e-14` for mean MCSE. All values
met the original `1e-8` absolute-plus-relative comparison tolerance.

The original serialized reader rebuilt every row from raw traces and reproduced
all source partitions. Closure checked 50,449 source/artifact bindings and both
original completion journals. The [execution checkpoint](../metadata/independent_baliphy_scalar_execution_checkpoint_20261002.json)
checks the three successful version-2 invocations and preserves the three
failed version-1 invocations. It rechecks compact archive/receipt hashes without
claiming to repeat the closure's complete large-data hash passes.

No quartet passes every scalar or candidate-length mixing screen at either
cutoff. This establishes reproducible marginal numerical calculations, while
adequate mixing, model adequacy and an accepted ancestral posterior remain open.

## Preserved percentile-boundary failure

The first full comparison stopped after 44 quartet checkpoints on one tail-ESS
discrepancy. Its [failed producer and dependency journals](../metadata/independent_baliphy_scalar_v1_failure_20261002.json),
code, plan and partial outputs are preserved. NumPy's interpolated percentile
and the locked report's weighted type-7 percentile rounded a tied 5% boundary
one representable value apart: 151 versus 150 indicator memberships. Direct-lag
and FFT ESS agree for identical indicators; the discrepancy arose before
autocorrelation estimation. It was not evidence of converged sampling.

Version 2 independently computes the weighted interpolation order described by
[SciPy's type-7 formula](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.mstats.mquantiles.html),
preserving the original comparison tolerances and diagnostic thresholds.
The exact source trace was added to the numerical fixtures and reproduced from
all four hashed raw logs with byte-identical array output. This is a software
regression case; the new production version still covers all 405 quartets.

## Reproduction and outstanding work

Run the producer and serialized reader with `--plan` and new explicit output/
launch identities. The [three-stage launch inventory](../metadata/independent_baliphy_scalar_launches_20261002_v2.json)
preserves the original invocations. For software checks, run
`check_independent_ancestral_scalar_diagnostics_v2.py --output NEW.json` in the
locked ancestral-diagnostic environment and
`check_independent_baliphy_scalar_readback_v2.py --output NEW.json` in the main
analysis environment.

The [main-environment reproduction specification](../environments/independent-ancestral-verification-20261002.yml)
records Python 3.10.13, NumPy 2.2.6, SciPy 1.15.3 and psutil 7.2.2. The
[observed backend snapshot](../metadata/independent_ancestral_verification_environment_20261002.json)
was recorded after closure; it is not a retroactive launch-time attestation and
does not alter frozen launch pins. Oracle tests separately require the existing
locked ArviZ environment.

Recreate the external regression array with
`materialize_scalar_quantile_boundary_regression.py --failure` using the
preserved failure record and a new `--output` directory. The fixture loader's
expected array location is recorded in the numerical-check script and full plan;
the plan also pins the array and its raw-source locator JSON.

Full production scalar numerical agreement/readback/provenance closure passed.
[Separate categorical software checks and a full source inventory](independent-ancestral-categories-20261002.md)
are now available, and full production categorical replay is running. Native alignment-parser checks, adequate
sampling horizons and joint posterior/model/root qualification remain open.
The two unresolved native allocation failures remain failures. Ancestral
structure prediction remains subject to GPU authorization; GPU inference is
paused. Aim 8 and all other biological aims remain unfinished.
