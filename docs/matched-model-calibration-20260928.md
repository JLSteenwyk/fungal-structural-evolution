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
