# Execution access restored: September 28 checkpoint

Shell access worked again at the September 28 11:38 EDT check. Background
analyses had continued during the access failure; no inference jobs were
restarted to recover observation access.

The matched comparative-model run completed all 144,040 unique fits. The
original output audit retains 102,975 optimization-review flags. The separate
analytic-gradient assessment and its readback completed for every fit,
leaving **658 unique fits** for review: 257 gradient-threshold cases, 395
optimizer-termination cases and six full-face start-disagreement cases.
The distinction matters: the original finite-difference review flags are
preserved, while the analytic check supplies additional numerical evidence.
Neither is a global-optimum proof or a calibrated biological test.

The full setting export contains 414,720 setting/tree rows, including 2,099
rows linked to those 658 unresolved fits. Refinement and inferential
calibration remain required. The separate 432,120-fit ordinary-ML polynomial
grid remains active.

Both FastML variants and the independent probability replay are terminal:
306 nonempty fits (153 per variant), plus six explicitly empty inputs.
Each variant has 8,058,340 independently replayed probability rows. The
maximum absolute fitted-parameter likelihood discrepancy is 0.473947975 for
precision-only and 0.003451690 for precision plus cache refresh. Respectively
77 and 74 fits exceed an absolute discrepancy of 1e-5; this is a descriptive
count, not a newly adopted acceptance threshold. Both variants have maximum
probability discrepancy 6.80132e-5. The largest raw probability-boundary
excursion is 2.66454e-15, retained without clipping. Remaining discrepancies,
optimization and model adequacy still require review.

Seven completed services were checked for inactive/success/exit 0. Receipt
bindings and 9,819 artifact files were verified, including nested FastML
attempt artifacts. The check does not rerun all optimizations or independently
recompute analytic derivatives. Evidence is in
`metadata/restored_access_terminal_checkpoint_20260928.json` and can be
reproduced with `scripts/record_restored_access_checkpoint.py --output NEW.json`.

The ancestral producer had 161 saved-sample audits and 161 successful state
extractions when checked, with 16 inference workers showing CPU progress.
This is not a convergence assessment. GPU structure prediction remains under
the prior pause; no GPU launch was made during recovery.

The interrupted observation of the categorical performance benchmark was
also recovered. All 120 synthetic patterns completed; their source hashes
passed verification. Full-grid diagnostic scheduling remains to be planned
from these measurements, not from an assumed benchmark failure.
