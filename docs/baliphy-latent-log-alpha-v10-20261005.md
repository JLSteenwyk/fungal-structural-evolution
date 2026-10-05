# Retaining latent log-alpha in future ancestral runs

The V10 logger retains the actual sampled log-alpha, its exponentiated alpha,
four rate categories, the Laplace prior parameters and its latent log density.
It passes all **405 reversible generated-source checks**, three paired native
prior controls and independent reconstruction of all **120 diagnostic rows**.
Every one of the **18 original scientific output files** is byte identical
between the paired V7 and V10 controls. Original tool, native and complete
invocation-journal evidence closes for both qualification and readback.
[Completed evidence](../metadata/baliphy_log_alpha_v10_completed_20261005_v1.json).

This is software qualification for the existing inference workflow. It does
not establish the cause of the historical infinite-alpha observations or an
adequate ancestral posterior. The [complete 24-role full-input V10 comparison](baliphy-latent-log-alpha-full-grid-20261005.md)
has now launched after construction and reader-software qualification. A new
full MCMC horizon has not launched. Existing reviewed ancestral arrays remain
excluded.

## Preserving the random operation and scientific outputs

The installed exponential distribution transform uses `exp <$> sample dist`.
Its `Random` functor implements `fmap` through a nested `RanBind`. V10 keeps
that nested operation while retaining both values:

```haskell
(alpha_2, projectLatentLogAlphaV10) <-
  (\draw -> (exp draw, draw)) <$> sample (laplace mu scale)
```

The three existing priors retain their original parameters: package `(6,2)`,
centered `(0,1)` and broad `(0,2)`. The existing scalar logging action writes
the additional `C1.P1.log-alpha-samples.jsonl`; the original scientific schema,
likelihood, seeds, priors and native caps are preserved. The diagnostic stores
explicit nonfinite tags and numeric quality counts rather than silently
replacing an overflow value. It does not derive a latent value from infinity.

The paired controls use 20 iterations and the same original seeds for each
prior. All six existing files agree byte for byte per prior: scalar TSV,
scalar JSON, column map, runtime tree, ancestral FASTAs and joint site-property
samples. All 63 diagnostic rows agree with their corresponding original
scalar/context records. The native controls retain their 12 GiB address-space,
300 CPU/wall-second and 64 MiB per-file limits.
[Source transformation](../scripts/baliphy_log_alpha_logger_v10.py) and
[qualification](../scripts/check_baliphy_log_alpha_logger_v10.py).

## Retained unsuccessful candidates

V8 and V9 used a direct latent draw instead of preserving the nested functor
operation. Their first paired broad-prior native control exited zero, but the
original scalar trajectories changed from iteration one. Both qualification
wrappers therefore exited one. The remaining two prior controls and the
deterministic probe were not run for those versions. No success receipt is
claimed for either candidate; all original sources, changed outputs, actual
tool failures and whole journals are retained.

V8 used a separate diagnostic action; V9 used the existing scalar action.
Their scientific and diagnostic outputs agree with one another for this
control. That observation alone does not establish a general causal claim
about random operation ordering. V10's retained nested operation passes the
unchanged-output test under all three tested priors.
[Complete retained failure evidence](../metadata/baliphy_log_alpha_noninterference_failures_20261005_v1.json).

## Independent diagnostic reconstruction

The separate reader imports none of the transformation, qualification or
original scalar-reader logic. It uses strict standard JSON, rejects duplicate
keys and literal nonfinite JSON numbers, and computes exponentiation and
Laplace log density with 90-digit Decimal arithmetic. It reconstructs all
63 native rows, the 18 paired file comparisons and **57 deterministic cases**:
19 latent values under each of the three priors. Fifteen deterministic rows
retain positive infinity in the derived alpha with finite latent values.

Ten altered semantic records are rejected, covering latent/derived values,
prior location/scale, density, rate normalization, leaf counts, Boolean values,
missing overflow tags and the wrong overflow state. Existing numeric
tolerances are unchanged. Rate normalization and the unit-rate large-shape
observations do not newly prove general Gamma-category accuracy. The complete
405-source reversibility test is producer evidence; the separate reader does
not independently compile or reverse all 405 programs.
[Independent reader](../scripts/readback_baliphy_log_alpha_v10.py) and
[its result](../metadata/baliphy_log_alpha_v10_readback_20261005_v1.json).

Qualification used two CPUs, 16 GiB RAM, no swap, 12 GiB address space and one
BLAS thread; the stage CPU/wall caps were 1,200/1,800 seconds. The independent
reader used the same CPU/RAM limits and 600/900-second stage caps, without new
MCMC runs. These are resource limits, not posterior run lengths or ETAs. No
GPU, paid infrastructure, original-job restart or installed-library edit was
used.

## Reproduction and remaining inference work

The large native outputs remain outside Git at
`data/software_audits/baliphy-log-alpha-v10-20261005-v1/`. Qualification and
reader resources, original tool payloads, source hashes, execution receipts
and complete journals are tracked in metadata. A reproduction needs the
pinned original native software and source fixtures, and fresh output paths.

The complete check against all 24 original full-input V7 comparison roles is now running,
preserving inputs, seeds, priors and caps and retaining every disposition.
Then resource-estimate full ancestral sampling and establish convergence and
adequate posterior uncertainty before biological admission. The installed
original formatter is not repaired. The historical covariance discrepancy,
accepted phylogenetic frameworks, structural accuracy calibration and all
eight evolutionary aims remain unresolved or incomplete.
