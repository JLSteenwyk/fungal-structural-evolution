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
