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
