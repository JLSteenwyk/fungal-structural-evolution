# Full-grid analytic stationarity assessment

Following the verified frozen-prefix diagnosis, a separate assessment is queued
for every one of the 144,040 production tree-fit dispositions. It waits for the
exact recorded output-integrity auditor to exit successfully and requires its
receipt to match the completed production receipt. The running optimizer,
original result files and original flags remain unchanged.

The assessment reconstructs each unique record input from the original
source partition and selection records. It rechecks both numeric and ordered
node/family/species identity fingerprints, recreates the scaled fixed design
and tree-factor rows, and checks the analytic objective against the saved fit.
It then computes the independently fixture-validated analytic derivative in
log(1 + variance ratio) units, applies the original boundary projection rules,
and compares with the unchanged maximum absolute gradient threshold of 1e-3.

Outputs distinguish the original gradient flag, analytic gradient flag, and
the remaining optimizer checks: successful optimizer termination, agreement
between full-face starts and absence of upper-bound contact. Passing an analytic
score alone cannot clear other numerical problems. Original execution errors
and analytic cases that cannot be evaluated reliably are retained explicitly.
No subset is dropped and no original status is overwritten. A passing local
stationarity check is not proof of a global optimum, valid uncertainty or an
adequate biological model.

The implementation processes one source partition at a time rather than holding
all reconstructed model data in memory. Source and script pins are checked
before and after execution. The final disposition set must equal the complete
production set, and each output row records its immutable source-fit hash.
Changing source data, scripts or dependency identities causes an explicit
failure rather than an implicit replacement or restart.

Script: `scripts/assess_full_working_model_gradients.py`.
Plan: `metadata/full_working_model_analytic_stationarity_plan_20260927.json`.
Launch: `metadata/full_working_model_analytic_stationarity_launch_20260927.json`.
Output: `results/model_validation/full-working-model-analytic-stationarity-20260927-v1`.

The job uses one CPU, 16 GiB memory, no swap and at most 1 GiB planned output.
The compute allowance is 2–48 hours after the production output audit, not a
measured completion ETA. At launch it is waiting on that audit. Full-grid
results, their readback, any remaining optimizer refinement, model adequacy
and inferential calibration remain pending. GPU prediction remains paused.

## Full output readback and refinement queue

A separate readback now waits on the analytic assessment's recorded process
identity and successful terminal service state. It verifies all 144,040 expected
recipe/tree combinations against the production manifest and immutable fit
bytes, then independently reconstructs the scalar boundary projection and
all convergence flags. Assessment counts and transitions must exactly match
the producer receipt. Duplicate, extra or missing dispositions fail the check.

Every fit error, unresolved derivative, failed analytic threshold, unsuccessful
optimizer termination, upper-bound contact or disagreement among full-face
starts is retained in `remaining_review_cases.jsonl`, with source hashes and
explicit reasons. This provides the complete next-stage refinement inventory;
it does not change original fits or clear other flags when a gradient passes.
The numerical derivative itself is not recomputed by this output readback.

Five valid boundary/error cases and ten altered cases exercised exact-threshold
handling, opposing boundary signs, nonfinite/malformed derivatives, false
clearance and preservation of original errors. All expected rejections passed:
`metadata/analytic_stationarity_readback_checks_20260927.json`.
These fixtures do not establish that the full production output has passed.

Run `scripts/readback_full_analytic_stationarity.py --plan
metadata/full_analytic_stationarity_readback_plan_20260927.json`. The queued
service identity is recorded in
`metadata/full_analytic_stationarity_readback_launch_20260927.json`.
Output goes to
`results/model_validation/full-working-model-analytic-stationarity-readback-20260927-v1`.
Resources are one CPU, 4 GiB memory, no swap, and 1 GiB planned output;
the 0.1–6 hour planning allowance excludes waiting for the assessment.
The full readback and remaining-review inventory are pending. Further
optimization, model adequacy and inferential calibration remain open.
