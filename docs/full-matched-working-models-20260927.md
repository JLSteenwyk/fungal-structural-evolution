# Complete matched working-model fitting run

The restartable production controller is queued for all 28,808 exact unique
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
