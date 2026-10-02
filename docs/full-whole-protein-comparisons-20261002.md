# Complete whole-protein comparison integration — October 2, 2026

The complete model-comparison workflow is implemented and queued for **75,070
unique numerical inputs × five tree alternatives = 375,350 fits**, with
**829,440 comparisons per tree = 4,147,200 comparisons**. It advances the
sequence–structure coupling and duplication working models (aims 2 and 4).
Production results have not completed. All eight biological aims remain open.

## Source cohort and selection

This version uses the existing frozen whole-protein measurement cohort. It does
not transfer these fits to the expanded structural atlas or new background
measurements. The older queued export that uses five explicit recovery
candidates remains a separate immutable version.

The prerequisite is the [complete optimization handoff](full-whole-protein-optimization-closure-20261002.md):
all original fits, input hashes, five-fit audit shards, every dynamic flagged
follow-up and its independent numerical reader, with four original
invocation-linked completion/resource journals and complete source hashes.
The new source loader requires that closure, the full Cartesian grid, exact
follow-up keys equal to all original flags, verified comparison links and the
resolved zero-reference-support registry.

| Original state / follow-up | Selected fit | Effective state |
|---|---|---|
| Original numerical pass | Original production fit | Original pass retained |
| Original fit error | Original error; null objective and parameters | Error retained |
| Original flag; passing closed refinement | Follow-up fitted result | Numerical refinement pass |
| Original flag; passing closed recovery | Follow-up fitted result | Numerical refinement pass |
| Original flag; unresolved closed refinement/recovery | Follow-up fitted result | Full-grid follow-up review required |
| Original flag; follow-up error | Original numerically audited flagged fit | Original optimization flag retained |

The last row preserves the original candidate, objective and parameters; the
failed follow-up's path, hash and status remain explicit. A numeric audit of a
candidate does not establish optimizer convergence or global optimality.

## Export and independent readback

`fit_summaries.jsonl` records original and effective states, original and
selected source paths/hashes, follow-up provenance and resolved support. The
normalized `selected_fit_parameters.jsonl.gz` contains one record per fit:
original column/specification identity, covariate scales, scaled and raw-unit
coefficients and conditional GLS coefficient covariance, selected log1p
variance ratios, profiled residual scale and ordinary-ML denominator. Original
errors retain null parameter fields. Every parameter record retains
`scientific_eligibility=false`.

The independent reader implements source selection separately from the new
producer. It reconstructs every source summary and parameter record, checks
dimensions and finite values, and independently verifies coefficient-unit and
covariance-unit transformations. JSON comparisons distinguish booleans,
integers and floating-point fields. The producer and reader share closed-source
I/O, the support registry and immutable prior likelihood comparison arithmetic;
the reader uses the prior separate arithmetic/orientation/decision verifier.
This is source-selection and numerical readback, not a second inference engine.

All comparison links are retained in their original order. Different
observations have no direct likelihood contrast; missing fits, unresolved
optimization, zero-reference-support failures and nested/identical-design
likelihood violations remain explicit. Negative likelihood gains are retained.
No p-values, chi-square calibration, model-selection decision, biological
acceleration or causal duplication effect is produced.

Conditional coefficient covariance describes the fitted Gaussian working model
at the selected variance parameters. It does not include variance-estimation,
tree, reconciliation, sampling, alignment, prediction or model uncertainty and
is not a calibrated confidence interval. Those controls and inference remain
required.

Parameter gzip headers have fixed time and no filename so regenerated
parameter hashes are identical on recovery. Each committed tree checkpoint
binds the frozen plan, complete scalar summary, parameter file and compressed
comparisons. Completed exports refuse restart. Full readback rechecks every
row, source hash, exact denominator, status count and tree checkpoint. Final
closure additionally requires both original producer/reader completion and
resource journals and all artifact/source hashes.

## Validation and resources

The [software validation locator](../metadata/full_whole_protein_comparisons_fixture_validation_20261002.json)
binds a ten-input/five-tree fixture: 50 summaries, 50 parameter records, 20
follow-ups and 225 comparisons across all five relationship types. All original
and follow-up pass/review/error, refinement and recovery paths passed. A run
interrupted after its first committed tree reused that checkpoint with an
identical parameter-file hash, and a completed restart was refused. All 26
invalid raw-source or rehashed export cases were rejected, including incorrect
coefficient/covariance transforms, dimensions, scales, promoted errors/reviews,
source identities, changed ratios, manufactured error coefficients, type
substitutions, missing rows, wrong arithmetic/nesting and incomplete follow-up
census. These are software contracts; production closure I/O is replaced only
in the fixture. Ten current original passed fits separately confirmed parameter
format compatibility; this is not full-grid verification.

The [prelaunch resource estimate](../metadata/full_whole_protein_comparisons_resources_20261002.json)
specifies two CPU, 32 GiB RAM, no swap, one BLAS thread, 64 GiB output and 100
GiB free reserve. A planning allowance of 1–48 active hours for each producer
or reader is uncalibrated and excludes prerequisite waits; it is not an ETA.
There are no GPU or paid resources. Existing jobs were neither changed nor
restarted; GPU prediction remains paused.

The [80-pin frozen plan](../metadata/full_whole_protein_comparisons_plan_20261002.json)
and [original launches](../metadata/full_whole_protein_comparisons_launches_20261002.json)
queue producer PID 1810771, reader 1810866 and closure 1810908 behind original
optimization closure PID 246680. The [fresh runtime checkpoint](../metadata/project_runtime_checkpoint_20261002_v19.json)
verified 51 pipeline/six original handles, 35 terminal successes, ten preserved
historical pipeline failures and 467,128 distinct closed bindings. The
[broader inventory](../metadata/project_live_launch_inventory_20261002_v2.json)
rechecked all 66 original live handles, including nine legacy-schema wrappers.
At that checkpoint, background alignment was 293,184/298,848; shortly afterward,
native sequence alignment was 195,712/324,672 and the original whole-protein fit
manifest had 188,365/375,350 dispositions. These are observed production counts,
not completion forecasts or biological acceptance.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 /home/bizon/anaconda3/bin/python \
  scripts/check_full_whole_protein_comparisons.py \
  --output /tmp/full-whole-protein-comparison-software-recheck.json

# After complete optimization closure, use the frozen plan. Do not duplicate
# the original queued jobs or overwrite their output directory.
/home/bizon/anaconda3/bin/python scripts/export_full_whole_protein_comparisons.py \
  --plan metadata/full_whole_protein_comparisons_plan_20261002.json
/home/bizon/anaconda3/bin/python scripts/readback_full_whole_protein_comparisons.py \
  --plan metadata/full_whole_protein_comparisons_plan_20261002.json \
  --output results/structural_comparisons/full-whole-protein-comparisons-20261002-v1/readback.json
```

Large outputs and complete proof maps remain outside Git. Compact plans,
validation/launch locators and reproducible scripts are versioned. Final
completion evidence is expected at
`metadata/full_whole_protein_comparisons_completed_20261002.json`; its absence
means this production stage has not passed final closure.
