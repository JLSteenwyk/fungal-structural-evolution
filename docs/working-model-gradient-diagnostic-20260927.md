# Separate diagnostic for working-model gradient flags

The first frozen 320-fit production prefix contained 180 projected-gradient
review flags, despite successful optimizer termination, agreement of all
full-face starts and no upper-bound contacts. The running fits and their review
criteria remain unchanged. A separate diagnostic examines every fit in that
preexisting prefix, including the 140 whose saved gradient check passed.

`scripts/matched_reml_gradient.py` differentiates the grouped cross-products,
nested family correction, low-rank species update, log determinants and profiled
residual scale analytically. It returns derivatives in variance-ratio units and
in the optimizer's log(1 + ratio) units. At a covariance V and residual projector
P, the independent reference identity for component derivative K is one half
of trace(P K) minus (n minus fixed-effect rank) times y-transpose P K P y divided
by y-transpose P y. `scripts/check_matched_reml_gradient.py` uses this dense
identity, rather than differentiating the cached expressions, to check all 243
component derivatives in 81 zero/moderate/large-ratio and factor configurations.
Maximum absolute ratio-gradient discrepancy was 7.03e-10. The analytic diagnostic
rejects cancellation-sensitive cross-products rather than silently reporting
an unreliable score.

`scripts/diagnose_working_model_gradients.py` reconstructs every input in the
frozen prefix from checksum-bound source recipes and verifies both numerical
and ordered-identity hashes. It compares the analytic objective to the saved
fit, computes finite differences with steps 1e-3, 1e-4, 1e-5 and 1e-6 in the
same parameterization, and checks that the 1e-4 result reproduces the saved
gradient. Analytic gradients are projected using the original boundary rules
and compared to the same 1e-3 diagnostic threshold, without changing the
original output status.

Derivative fixture proof: `metadata/matched_reml_gradient_checks_20260927.json`.
Diagnostic launch: `metadata/working_model_gradient_diagnostic_launch_20260927.json`.
Output: `results/model_validation/working-model-gradient-diagnostic-20260927-v1`.
The job uses one CPU, 16 GiB memory and no swap, with 0.1–2 hours allowed.

The full 320-fit diagnostic completed successfully, and all output rows, source
pins and terminal state were checked. All 140 original gradient passes remain
below the same threshold analytically. Of 180 originally flagged fits, 179 are
below it analytically and one remains above (maximum component magnitude
0.00109108 versus threshold 0.001). The median maximum component discrepancy
from the analytic score decreases from 0.02050 at finite-difference step 1e-4
to 0.0003207 at 1e-5 and 0.00001309 at 1e-6. Boundary-truncated finite differences
can retain larger errors, so simply choosing one smaller step is not validated
as a universal fix. The analytic and saved objectives agree within 1.82e-12.
Completion proof: `metadata/working_model_gradient_diagnostic_completed_20260927.json`.

This is a numerical investigation
of recorded flags, not a biological pilot, full-grid reclassification, proof of
global optimization, uncertainty calibration or permission to discard review
cases. Evidence of finite-difference error would motivate a separately verified
full-grid analytic assessment; it would not erase existing audit history.
