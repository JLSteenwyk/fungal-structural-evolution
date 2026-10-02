# Full expanded uniform covariance qualification

The complete stage covers all **622,080 original model-setting records**, both
signed/unsigned loading modes and all five working species-tree alternatives:
**6,220,800 setting/mode/tree links**. It is queued behind the original design-v2
provenance closure. No production result is accepted until the independent
reader and final full source/artifact/original-journal closure finish.

The original design producer reports 4,340 exact cohorts, implying 130,200
distinct fixed designs and 1,302,000 design/mode/tree audits. These counts are
conditional until complete design readback and closure. Exact reuse avoids
repeating a covariance calculation for identical source data/designs; every
original policy, scenario, mask, screen, outcome, degree and physical order
remains in the exported setting links. Empty and excluded records remain.

## Explicit uniform variance definitions

Every cohort directly checks the target-node incidence: one unit entry per
row in distinct columns. Its covariance is identity, so target-node variance
is combined with uniform residual variance. Signed endpoint-family loadings
cancel exactly. Unsigned endpoint-family incidence is twice the one-per-case
family contrast-intercept incidence, giving a kernel four times as large.
Their combined variance is `variance_intercept + 4*variance_endpoint_family`.
These exact combinations preserve the original nonnegative covariance cone.

The seven represented kernels are uniform residual identity, background-node,
physical model-pair, gene, versioned model, family contrast-intercept and
species covariance. Original distinct genes and shared model coordinates
remain separate namespaces. Every audit records the folded definitions.
Additional dependencies are review states; no arbitrary basis is deleted.
Nonuniform residual weighting and other controls require separate definitions
and qualification before final model identities or inference.

## Numerical and independent checks

Each active full-rank fixed design is reconstructed from closed input
partitions and checked against its original design and input-ID hashes.
Raw covariance kernels and their REML error-contrast restrictions `H K H`,
where `H=I-QQ'`, receive separate Gram, resolution and rank audits. Source
rank/empty states bypass numerical qualification explicitly. Constant
responses retain their original exclusion even if a covariance basis passes.
Numerical failures remain review records with their original error details.

The producer computes exact component-local entity products once per cohort
and loading mode and reuses them across all fixed designs and trees. The
independent reader uses sparse latent overlaps `Z_a' Z_b`, pivoted QR with
extended-precision column normalization, and independent `gesvd` checks.
For entity kernels it evaluates `||(H Z_a)' (H Z_b)||_F^2` through latent
overlap and low-rank correction identities. Species contractions are
reassociated to avoid caching dense latent-by-species matrices for every tree.
No global dense covariance, residual projector or latent-overlap matrix is
allocated. SQLite enforces unique audits and verifies every original setting
link and disposition. A producer-qualified basis must also pass independent
qualification; conservative review differences remain counted.

[Primitive validation](../metadata/covariance_basis_context_validation_20261002_v2.json)
passed 48 dense/block/latent comparisons, including signed loadings, empty
global columns, zero-rank species factors and extreme design units. Maximum
relative projected-Gram discrepancy was 3.12e-16. Seven invalid inputs were
rejected. [Full software contracts](../metadata/full_uniform_covariance_qualification_validation_20261002_v3.json)
passed all 720 synthetic settings, 1,800 design audits and 7,200 links, the
independent seven-basis qualified path, interrupted replay, completed-restart
refusal and all 17 rehashed altered-export rejections. Earlier unlaunched
validation receipts remain preserved; version 3 binds the final frozen
implementation. These checks are software fixtures, not a biological pilot.

## Resources and reproduction

Before launch recorded **two CPU cores, 32 GiB RAM, no swap and one BLAS
thread**. A 64 GiB output/scratch budget and 164 GiB minimum free-space check
cover exports, SQLite and complete provenance receipts. The worst pre-sharing
record budget is 53.4 GiB uncompressed; compression reduces actual output.
The largest five-kernel component stack is 332,698,240 bytes. Complete source
proof dictionaries may dominate memory. The 1–168 hour planning range per
stage is uncalibrated; full production timings will be reported as they arrive.
This is neither an optimization nor a project ETA. No GPU or paid resources.

The [production plan](../metadata/full_uniform_covariance_qualification_plan_20261002.json),
[resource record](../metadata/full_uniform_covariance_qualification_resources_20261002.json),
[dictionary](../metadata/full_uniform_covariance_qualification_dictionary_20261002.tsv)
and [environment](../environments/full-uniform-covariance-qualification-20261002.yml)
record exact identities, source closures, pins and execution requirements.

```bash
/home/bizon/anaconda3/bin/python scripts/launch_full_covariance_qualification.py \
  --plan metadata/full_uniform_covariance_qualification_plan_20261002.json
```

This version is queued and refuses duplicate launch/output identities.
Reproduction requires explicit new identities without overwriting original
jobs, plans or evidence.

## Scientific limits and next work

Qualification establishes numerical structural identifiability for a supplied
uniform working covariance. It is not a variance estimate, accepted ancestry
model, adequacy test, calibrated uncertainty or evolutionary association.
Expanded optimization, weight/control variants, uncertainty and multiple
testing, prediction-source/domain/PAE/ascertainment controls and accepted
species/reconciled gene phylogenies remain required. All eight aims remain
incomplete. Structural GPU prediction remains paused.
