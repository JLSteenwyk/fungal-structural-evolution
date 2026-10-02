# Complete provenance for whole-protein optimization follow-up

A final handoff is implemented and queued behind the original full numerical
follow-up reader. It covers **75,070 exact inputs across five tree alternatives,
375,350 original ML fit dispositions, all 75,070 five-fit audit shards and every
flagged follow-up**. Full production and numerical audit remain pending. This
supports the sequence–structure/duplication working models in the
[whole-protein workflow](whole-protein-workflow-20260929.md).

The original fit, audit, follow-up and reader scripts/plans remain immutable.
The initial five-candidate comparison export remains a separate frozen version.
Integration of the complete follow-up into a new full comparison export is still
required; this provenance handoff does not perform that integration or establish
calibrated biological effects.

## Exhaustive source and result checks

The handoff requires the complete 75,070-by-five Cartesian grid in both original
fit and audit manifests. Numerical passes, review flags and original fit errors
remain explicit. Every input file and original fit file is bound by its recorded
checksum. Each five-fit audit shard must have the exact original input, five
fit checksums and fit/audit plan identities; its result digest and all five
results must match the full audit manifest. Missing, duplicate or inconsistent
records cannot be represented as completed numerical evidence.

The follow-up scope manifest must equal the complete original disposition
census. Follow-up and independent-reader keys must equal the entire original
flag set. Each case retains original fit/input identities and checksums, its
result digest, source checksum and exact reader status. Passed, still-under-review
and failed follow-ups remain distinct. Original errors and follow-up errors
retain unverified numerical status; hashing an error does not qualify it as a fit.

The archive binds original input/source maps, all receipts/manifests, every
original fit/audit shard, every flagged case and all relevant pins. All four
original fit/audit/follow-up/reader processes must be absent, with exact
PID/command-linked original invocation completion/resource journals. A collected
unit's success status alone is insufficient. Complete numerical likelihood,
coefficient/covariance and gradient replay is performed by the existing full
audits; this final stage binds their complete evidence and checks linkage.

## Software checks and resources

The software contracts cover three inputs × five trees, all 15 original audit
states, three flags with passed/review/error follow-ups, and one retained
original fit error. All 11 corrupted evidence cases were rejected: missing and
duplicate audit results, changed input/fit hashes, false numerical qualification,
missing/duplicate reader records, changed reader source hashes, promoted
follow-up errors, dropped scope records and rehashed false original provenance.
These are implementation checks, not a pilot or full production result.

Resource estimates preceded launch: two CPU cores, 32 GiB RAM, no swap, 4 GiB
archive output and 100 GiB free-disk reserve. The minimum direct inventory has
525,490 input/fit/audit bindings; if every fit were flagged, it would have
900,840 direct bindings before inherited source maps and journals. A 1 KiB
per-entry planning scenario is about 0.86 GiB for that latter inventory.
The prior materialized-input resource estimate is about 62.75 GB. Arrays/fits
are hashed without copying or reoptimizing them. The 1–24-hour active planning
range is uncalibrated and excludes dependency waits; it is not an ETA or upper
bound. No GPUs, paid infrastructure or existing-job reconfiguration are used.

- [Frozen plan](../metadata/full_whole_protein_optimization_closure_plan_20261002.json)
- [Resource estimates](../metadata/full_whole_protein_optimization_closure_resources_20261002.json)
- [Passed software contracts](../metadata/full_whole_protein_optimization_closure_fixture_validation_20261002.json)
- [Original queued handoff identity](../metadata/full_whole_protein_optimization_closure_launch_20261002.json)
- [Full source/result handoff](../scripts/close_full_whole_protein_optimization.py)
- [Software checks](../scripts/check_full_whole_protein_optimization_closure.py)
- [Launcher](../scripts/launch_full_whole_protein_optimization_closure.py)

Large archives remain outside Git in
`results/model_validation/full-whole-protein-optimization-closure-20261002-v1/`.
The future completion locator is
`metadata/full_whole_protein_optimization_completed_20261002.json`.

Reproduce the software contracts with a new path:

```bash
python scripts/check_full_whole_protein_optimization_closure.py --output results/software-checks/whole-protein-closure-new-run/receipt.json
```

Fresh production requires new plan/output/unit identities and successful
completion of the original dependencies. Run
`python scripts/close_full_whole_protein_optimization.py --plan metadata/new_whole_protein_closure_plan.json`
only with those new pinned identities. Existing live jobs and completed evidence
must not be duplicated or overwritten.

## Corrected live-job inventory

A scan of every `metadata/*launch*.json` carrying PID/create/CMD fields, combined
with the previous original handles, verified **63 distinct live project jobs**.
Thirteen older jobs were absent from the recent 49-handle checkpoint; adding
them and this new handoff accounts for the increase. Four older whole-protein
follow-up/comparison wrappers have the current plan-bearing schema and were
added to the standard collector, now covering 48 pipeline plus six original
jobs. Nine older phylogenetic, polynomial and ancestral wrappers retain legacy
schemas and are recorded separately, with their original identity and available
command-plan/script hashes checked. No legacy launch record was rewritten.

- [Fresh 54-handle/full-hash checkpoint](../metadata/project_runtime_checkpoint_20261002_v18.json)
- [Complete 63-handle original live inventory](../metadata/project_live_launch_inventory_20261002.json)

These are live observations, not terminal scientific results. Numerical
stationarity, source integrity and covariate support do not establish global
optimality, identifiability, accepted phylogeny, model adequacy, calibrated
likelihood tests or biological inference. Full versioned comparison integration
and calibration remain required. All eight scientific aims remain incomplete;
GPU prediction remains paused.
