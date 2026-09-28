# Analytic refinement of remaining working-model fits

The completed 144,040-fit working-model grid has 658 unique tree fits that still
require numerical review: 257 fail the analytic projected-gradient threshold,
395 have optimizer-termination flags, and six have disagreement across full-face
starts. These represent 542 distinct record inputs. Original outputs and their
statuses remain preserved.

A separate refinement now runs every flagged fit. It reconstructs the exact
record order, values and species factors from the original pinned sources. It
uses the independently checked analytic REML derivative, all eight component
boundary faces and the original three nonempty-face starts. It also optimizes
from the original parameters and retains the original point as a candidate:
24 candidates per fit. A failed candidate with a lower objective remains visible.
Bounds and numerical acceptance thresholds are unchanged. The iteration limit
increases from 300 to 1,000; every candidate is checked against the separate
residual-based likelihood evaluator.

Validation used three synthetic covariance scenarios, checked all 72 candidate
likelihoods against dense algebra, and compared the best solutions with three
independent dense Powell searches per scenario. Maximum best-objective
disagreement was 2.05e-8. Four invalid parameter inputs were rejected. This is
implementation validation within the full workflow, not a reduced sampling study.

The production run uses one CPU core, at most 16 GiB RAM, no swap and no GPU.
The 658 original fits required 2,219.84 summed wall seconds; a conservative
1–8 hour allowance includes refinement and repeated source-table joins. This
is a planning range, not a measured completion forecast. Outputs are budgeted
at 1 GiB. No paid resources are involved.

Reproduce with:

```bash
OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python \
  scripts/run_matched_reml_analytic_refinement.py \
  --plan metadata/matched_reml_analytic_refinement_plan_20260928.json
```

The launch identity is recorded in
`metadata/matched_reml_analytic_refinement_launch_20260928.json`; the service is
`fungal-matched-reml-analytic-refinement-20260928.service`. Per-fit results are
restartable and bound to the plan and original source fit. Output location:
`results/structural_comparisons/matched-reml-analytic-refinement-20260928-v1`.
The fixture receipt is versioned in
`metadata/matched_reml_analytic_refinement_checks_20260928.json`.

At this checkpoint the production run is active, not complete. Next requirements
are terminal-state verification, independent output readback, review of any
remaining flags and explicit integration with the original setting map.
Numerical success does not establish component identifiability, calibrated
uncertainty, biological model adequacy or causal duplication effects.


The full output-readback service is now queued behind the refinement producer:
`fungal-matched-refinement-readback-20260928.service`. It revalidates the
producer's PID, creation time and command while waiting, then requires terminal
success before checking all 658 source-bound dispositions. It independently
checks all 24-candidate enumerations, fixed bounds, best-candidate selection,
projected-gradient arithmetic, status flags, coefficient-unit conversion and
explicit remaining-review records. This is not another derivative calculation
or model refit.

Before launch, the readback passed all 152 currently completed production
records and rejected altered candidate counts, gradient flags and raw-unit
coefficients. That frozen subset is documented in
`metadata/matched_refinement_readback_checks_20260928.json`; it is not a claim
about the remaining cases. The complete readback uses
`metadata/matched_refinement_readback_plan_20260928.json` and will write to
`results/model_validation/matched-reml-refinement-readback-20260928-v1`.
It allows one CPU core and 8 GiB RAM, with no GPU or paid resources.

Integrating qualified refinements with the complete original setting map,
model adequacy and inferential calibration remain subsequent requirements.


## Full estimate integration queued

`scripts/integrate_matched_refinements.py` waits for the identified refinement
readback service to reach inactive/success/exit0, then checks its receipt and
all bound artifacts before integrating all 658 targeted dispositions. The
new export preserves every original column in the 144,040 unique-fit and
414,720 full-setting rows. Separate `selected_*` columns carry the refined
coefficients, conditional intercept variance, variance ratios and objective;
source hashes, selection labels and review flags remain explicit. All 2,099
expanded rows linked to targeted fits are accounted for. An unsuccessful
refinement retains its original estimates and a review flag. Omitted
covariates stay null rather than becoming zero.

The mapping checks passed on 282 frozen completed refinements, including
unchanged original columns across the full setting table. Duplicate fits,
changed source hashes and worsened objectives were rejected; explicit
refinement-error fallback was checked. Frozen checks live under
`results/model_validation/matched-refinement-overlay-checks-20260928-v2`.
The first check invocation completed its assertions but failed when writing
metadata with a string instead of a Path; its v1 artifacts are preserved,
and the corrected v2 run is the recorded test result.

The queued integration has one CPU, 8 GiB memory, no swap and a 1 GiB storage
allowance; estimated active work is 0.02–1 hour after the upstream wait.
Plan and live launch identity are in
`metadata/matched_refinement_integration_plan_20260928.json` and
`metadata/matched_refinement_integration_launch_20260928.json`. These describe
the original queued attempt; the completed recovery is recorded below.
Conditional covariance is not calibrated
uncertainty, and numerical refinement does not qualify a biological effect.

## Completed refinement and integration

The version 2 output auditor and integration jobs both reached
inactive/success/exit0. The completed export is
`results/structural_comparisons/refined-working-model-grid-export-20260928-v2`.
All 658 targeted fits use audited refinement candidates; the other 143,382
unique fits retain their original estimates. These map to 2,099 refined and
412,621 unchanged setting rows, respectively. No selected numerical-review
flags remain. Original estimates and flags are preserved in the export.

`metadata/matched_refinement_integration_completed_20260928_v2.json` records
the full integration readback and receipt hashes. The export receipt and both
Parquet artifact hashes were rechecked against those records when updating
this document. Earlier failed attempts remain available in the job history.
The queued-stage descriptions above are historical, not current completion
claims. Numerical acceptance does not establish global optimality, calibrated
uncertainty, covariance adequacy or a biological effect. Residual replay and
candidate interval analyses are tracked in
[matched-model calibration](matched-model-calibration-20260928.md).
