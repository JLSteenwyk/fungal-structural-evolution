# Complete matched working-model fitting run

The restartable production controller is running all 28,808 exact unique
record inputs and all five species-tree alternatives: 144,040 tree-fit
dispositions. The complete source map retains all 82,944 original settings
(414,720 setting/tree results). Reuse must pass the full independent inventory
audit before any model is fitted. This is the complete specified linear working
model grid, not a favorable subset or a biological pilot.

The response is the target-minus-background RMSD difference averaged equally
over eligible domains within each selected record. Four domain-averaged
covariate differences enter the linear fixed design. Zero constant columns are
explicitly omitted; nonzero constant columns fail rather than silently changing
the meaning of the intercept. Varying covariates are divided by their population
standard deviations without centering, preserving the intercept at zero
covariate differences. Coefficients and conditional covariance are also
exported in original units. These conditional quantities are not final
confidence intervals or causal duplication effects.

Each model includes positive residual variance, nested shared-background and
family-component terms, and the specified additive endpoint species covariance.
The cached optimizer retains all 22 face/start attempts, boundaries, gradients,
convergence messages and review flags. Every candidate likelihood is recomputed
with the separately checked direct evaluator. Best coefficients, profiled scale,
residual quadratic form and conditional covariance are also checked directly.
Numerical errors are saved with their tracebacks. More than ten percent errors,
after at least ten errors, stops the controller for review rather than allowing
a pervasive implementation problem to consume the whole grid. Optimization
review flags alone are retained as scientific diagnostics.

Every submitted input is ordered by target node and rechecked against the
audited numeric and identity fingerprints. Per-tree results are written
atomically and include source-plan and payload hashes. Restart requires the
same plan and audit binding, validates existing payloads, and reuses only
matching saved results. A process lock prevents simultaneous controllers in
the same output directory. The work queue is bounded to twice worker count.
The final disposition manifest retains each output checksum and status; final
completion requires all expected inputs and tree alternatives, regardless of
whether their status requires subsequent review.

`scripts/check_full_matched_model_worker.py` passed five synthetic tree cases
and 110 direct candidate checks, constant-zero removal, exact restart, and
rejection of a corrupted payload. These fixtures check execution and restart
semantics, not model adequacy or inference calibration.

Controller: `scripts/run_full_matched_models.py`.
Plan: `metadata/full_matched_working_model_plan_20260927.json`.
Launch: `metadata/full_matched_working_model_launch_20260927.json`.
Worker proof: `metadata/full_matched_model_worker_checks_20260927.json`.
Output: `results/structural_comparisons/full-matched-working-models-20260927-v1`.

Resources are eight one-thread workers under a total eight-CPU quota, 48 GiB
memory, no swap and a 10 GiB output allowance. The planning range is 12 hours
to two weeks, based on measured evaluator costs with unknown production
optimizer behavior; this is not a measured completion ETA. The controller
first waits on the exact recorded inventory-verifier process and requires its
successful authoritative terminal state and matching audit receipt. CPU fitting
is authorized; GPU prediction remains paused and no paid resources are used.

All output initially remains a working-model result pending full run audit,
optimizer-case review, covariance adequacy, nonlinear sequence sensitivity,
prediction uncertainty, resampling/calibration and multiplicity handling.
Completing this run does not complete the project or its duplication aim.

The source audit has passed and production fitting has begun. A frozen
complete-line prefix of the manifest records 320 dispositions: 140 passing
numerical optimizer checks and 180 requiring review, with no execution errors.
Every review flag in that prefix is a projected-gradient failure; all full-face
starts agreed, optimizer statuses were successful, and no upper variance bound
was reached. `metadata/full_working_model_early_diagnostics_20260927.json` binds
the exact prefix and includes examples. This ordered early subset is not
representative evidence for the final review rate or biological effects. The
flags remain unchanged pending a separate gradient-accuracy investigation.

`scripts/audit_full_matched_model_outputs.py` is queued behind the exact live
controller. After successful terminal completion, it checks all expected
dispositions, file/payload/source checksums, the complete face/start enumeration,
parameter bounds, review classification, scale/degrees of freedom and raw-unit
coefficient transformations. It verifies the accounting of producer numerical
readbacks without claiming to repeat those numerical fits independently.
Launch: `metadata/full_matched_working_model_output_audit_launch_20260927.json`.
Expected proof: `metadata/full_matched_working_model_output_audit_20260927.json`.
This audit uses one CPU, 8 GiB memory and no swap, with 0.1–12 hours allowed after
production finishes. Error and review dispositions remain explicit.

The complete [joint-covariate support check](joint-covariate-support-20260927.md)
has passed for all 28,808 unique inputs and all 82,944 settings. Independent
readback accepted 28,784 original certificates; a separately documented
nonnegative-weight construction supported the remaining 24 within the unchanged
1e-8 scaled tolerance. This clears the numerical reference-support diagnostic,
not the pending covariance, optimization, uncertainty or biological gates.

## Complete reporting grid queued

After terminal success of the analytic-stationarity readback,
`scripts/export_matched_model_grid.py` will export all 144,040 unique fits and
all 414,720 original setting/tree combinations. Every fit retains original
status, analytic assessment, explicit numerical-review reasons, source-file
hash, coefficients in original units, conditional intercept variance, profiled
scale, objective and variance ratios. Omitted covariates remain null, never
zero. Failed fits retain their inventory record counts and missing estimates.
All original and projected joint-support dispositions remain joined to settings.

The exporter requires the full source-receipt chain and hashes, recomputes review
flags through the separately checked scalar logic, verifies the complete
input-by-tree key set, checks record counts, and reads both Parquet exports back.
Every expanded field must be constant across settings reusing a unique fit.
Native coefficient, omitted-covariate and failed-fit retention checks passed;
these are export checks, not model validation or calibrated inference.

Plan/launch: `metadata/matched_model_grid_export_{plan,launch}_20260927.json`.
Output: `results/structural_comparisons/full-working-model-grid-export-20260927-v1`.
The queued stage has one CPU, 8 GiB RAM, no swap, a 2 GiB output allowance, and
0.1–8 hours planned after the upstream readback. No confidence intervals,
p-values or biological effects are inferred by the export. Full-grid execution
is pending; this creates the reporting input for subsequent model assessment.
