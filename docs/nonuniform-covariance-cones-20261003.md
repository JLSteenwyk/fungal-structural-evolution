# Full covariance-cone proof for a separate residual diagonal

All 4,340 cohorts and 8,680 signed/unsigned certificates have completed an
exact covariance-cone proof and independent full readback that keeps the
residual diagonal separate from target-protein identity covariance. Complete
source/artifact/two-original-journal closure binds 4,429 hashes. This removes
an algebraic obstacle to nonuniform residual sensitivity models. It does not
choose weights, qualify their numerical identifiability, calibrate a
confidence-to-variance relationship or fit weighted evolutionary effects.

## Why the uniform fold cannot be transferred

The original target-node kernel equals the identity matrix in every closed
cohort. With uniform residual variance, both terms use that same matrix and
their nonnegative variances can be combined. For a general positive diagonal
`D`, the residual kernel is `D` while the target kernel stays `I`. Combining
them would change the covariance model.

The new proof retains both terms. Exact closed operator identities still give
`gene = 0.5 × (I + background)` and `model = 0.5 × pair`. The pair kernel
equals `I + background` in 3,616 cohorts per mode, while 724 genuine exceptions
keep it. Signed family covariance is zero; unsigned family covariance is
four times the family-intercept covariance. These existing integer operator
certificates are reused unchanged, without numerical rank-based deletion.

| Certificate group | Certificates across both modes | Uniform retained terms | Separate positive diagonal retained terms |
| --- | ---: | ---: | ---: |
| Exact pair identity | 7,232 | 4 | 5 |
| Pair exceptions retained | 1,448 | 5 | 6 |

The retained terms are residual `D`, target `I`, background, optional
model-pair, family intercept and species. The forward map combines each
original variance using nonnegative dyadic coefficients. Its nonnegative
right inverse assigns retained coefficients to the corresponding original
terms and sets the redundant original coefficients to zero. Their product
equals the identity exactly, proving equality of covariance cones for any
finite strictly positive diagonal. The independent reader reconstructs each
original column image with Python `Fraction`, rather than repeating the
producer's coefficient updates.

The result also holds for constant diagonals, but the extra term is then
dependent. Nearly constant diagonals may cause numerical reviews. No
algebraic certificate asserts raw or REML identifiability. Weighted Grams,
whitening, roundoff envelopes and rank/conditioning diagnostics must be
computed anew once actual diagonals are defined. Uniform numerical envelopes
are not inherited by this proof.

## Qualification and execution

Software qualification used both pair identities and loading modes, with
five residual patterns: one, a different positive constant, alternating
heterogeneity, a broad positive ramp and near-uniform values. It passed:

- 610 forward/right-inverse dense covariance roundtrips across 20 cases;
- 120 ML/REML production/dense/component-spectral and row-permutation checks;
- 360 finite-difference score cells at interior and zero-variance boundaries;
- 56 altered certificate/map rejections and 20 invalid diagonal rejections.

Maximum dense covariance error was 1.42e-14, objective error 5.68e-14 and
gradient error 8.53e-14. Finite-difference error reached 6.56e-6 and passed the
unchanged combined relative/absolute score policy, 3e-5/3e-7. These are small
synthetic numerical checks, not a biological pilot or real-data fit.

The complete source/export/readback software fixture includes all five
declared nonempty cohorts and ten certificates. Twelve rehashed output
alterations and six rehashed source alterations were rejected; completed
producer/reader restarts were refused. Restored positive bytes and complete
positive source bindings passed. Its parent source/journals are explicitly
synthetic. The first attempt failed because the fixture parent directory was
missing. That original checker, native exit, stderr and journal remain
preserved. Fresh checker V2 creates only that directory; the mathematical
and production scripts stay unchanged.

Both original software waits exited zero, with exact wrapper PID/create/
command/invocation journal payloads matching their native execution receipts.
The numerical checker used 10.35 child CPU seconds with reported peak RSS
108,527,616 bytes; the full source checker used 8.95 seconds and 169,869,312
bytes. Manager memory accounting is distinct from child RSS.

The real stage then checked the entire closed parent archive and every
original certificate, preserving all 75,188 logical cases and 68,220,240
cohort-row occurrences. The independent reader checked every serialized
field, original certificate digest, exact coefficient map and right inverse.
Producer, reader and closer have actual original terminal-success evidence.
A fresh closure check rehashed all 4,429 source/artifact bindings and verified
all three original terminal handles. Original operators, certificates,
uniform numerical jobs and queued fitting/timing plans remain unchanged.

Evidence:

- [Numerical software qualification](../metadata/nonuniform_covariance_cone_software_validation_20261003_v1.json)
  and [actual original wait/journal proof](../metadata/nonuniform_covariance_cone_software_transport_20261003_v1.json).
- [Preserved directory failure](../metadata/full_nonuniform_covariance_cone_software_failure_20261003_v1.json),
  [fresh source/export qualification](../metadata/full_nonuniform_covariance_cone_software_validation_20261003_v2.json)
  and [actual original wait/journal proof](../metadata/full_nonuniform_covariance_cone_software_transport_20261003_v2.json).
- [Full immutable plan](../metadata/full_nonuniform_covariance_cone_plan_20261003_v1.json),
  [prelaunch resources](../metadata/full_nonuniform_covariance_cone_resources_20261003_v1.json)
  and [original launch handles](../metadata/full_nonuniform_covariance_cone_launches_20261003_v1.json).
- [Full completed proof](../metadata/full_nonuniform_covariance_cones_completed_20261003_v1.json)
  and [fresh complete provenance verification](../metadata/full_nonuniform_covariance_cone_completed_verification_20261003_v1.json).

## Resources, reproduction and next gates

The full metadata proof used two CPUs/16 GiB/no swap, one BLAS thread, native
12-GiB address-space, 3,600-CPU-second and 256-MiB-per-file limits. The
prelaunch record allowance was 16 KiB per certificate, 135.625 MiB total,
with 1 GiB output allowance and 64 GiB minimum free disk. The provisional
1–600-second stage range was uncalibrated. These are resource allocations
and record budgets, not measured process peaks or a total directory quota.
Existing authorized local CPU resources were used, with no new charge or
GPU work.

The versioned environment is the existing shared-entity fitting environment.
`check_nonuniform_covariance_cone.py` runs numerical qualification;
`check_full_nonuniform_covariance_cones_v2.py` runs source/export qualification.
Reproduce each under its saved wrapper/resource configuration using new output
and receipt paths. `run_full_nonuniform_covariance_cones.py` has separate
producer and `--reader` modes. The immutable plan and original recorded
commands bind the exact full source, resources and scripts. Restore the
pinned versions and source archives before reproducing on another machine;
do not restart a completed or failed original root.

Actual residual diagonals still need complete-case provenance and justified
normalization, weighting interpretation and missingness controls. Predictor
confidence scores are not calibrated error variances. A future sensitivity
policy must be identified as an assumption unless separately validated.
Every weighted design requires fresh raw/REML numerical qualification,
appropriate complete-scope timing, enforced resources and independent fits.
Inference/calibration, nonuniform controls, accepted phylogeny/reconciliation,
ancestral posterior adequacy, the full atlas and all eight biological aims
remain incomplete. GPU inference stays paused.
