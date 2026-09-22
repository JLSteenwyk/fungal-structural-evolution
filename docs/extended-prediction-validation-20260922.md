# Completed extended marker predictions and validation

The 769–1,024-residue prediction queue finished on September 22, 2026.
Its 4,363 candidate sequences are present in three disjoint output directories:
2,318 in the original interrupted run, 1,022 in the GPU 0 continuation, and
1,023 in the GPU 1 continuation. The original interruption receipt remains
unchanged; its remaining 2,045 sequences were completed by the two continuations.
The separate 513–768-residue cohort contains 5,510 audited, converted predictions.
These are marker cohorts, not predictions for the full fungal proteomes.

Independent validation and mmCIF conversion of all 4,363 new predictions
**passed** on September 22. The completed combined receipt is archived as
`metadata/extended_partitioned_handoff_20260922_completed.json`. The CPU-only controller is
`scripts/advance_partitioned_prediction_snapshot.py`, configured by
`metadata/extended_partitioned_handoff_20260922_config.json`. The service is
`fungal-extended-validation-20260922.service`; its observed process identity is
recorded in `metadata/extended_partitioned_handoff_20260922_launch.json`.

Before starting, the controller verified disjoint receipt identities and exact
coverage of the original FASTA. It runs the established artifact audit and
mmCIF conversion sequentially for each partition. Validation checks source
hashes, sequence and residue numbering, finite atomic coordinates, pLDDT and
PAE arrays, and lossless mmCIF roundtrips. Each partition keeps its original
prediction configuration and model provenance. The final gate also checks the
exact union of audited sequence IDs and the full taxon–marker link multiset,
including duplicate link multiplicities. The original partial output is
explicitly audited as a partial snapshot; only the combined union can establish
completion of the original queue. Six previously reused sequences remain
outside this 4,363-sequence prediction queue.

Resources are capped at two CPU equivalents and 16 GiB RAM, without swap or
GPU visibility. The planning allowance is 0.5–6 hours and at most 50 GiB of new
output, with a 100 GiB free-disk gate. These are planning estimates, not measured
completion times. At launch, approximately 12 TiB disk and 944 GiB RAM were
available. No paid resources were provisioned.

Logs and the eventual completion receipt are under
`results/predictions/extended-partitioned-handoff-20260922-v1/`. Each partition
has separate immutable audit and converted-structure directories listed in the
configuration. Failed or interrupted execution requires inspection and fresh
output paths; rerunning cannot overwrite an existing snapshot.

Successful handoff establishes artifact integrity and serialization, not
experimental structural accuracy. Residue mapping, PAE binding, native feature
extraction, matched predictor comparisons, and phylogenetic integration remain
separate downstream tasks. GPU prediction remains stopped after the authorized
run window.

## Whole-cohort residue integration

The 4,363 candidate sequence identities match 4,384 global marker-protein links
across 99 markers and 276 taxa. The subsequent controller verified the successful combined validation receipt
and began mapping. It combined the three audited inventories without copying coordinates or
changing model IDs, prediction configuration hashes, original receipt references
or confidence-array paths. The original interrupted partition remains explicitly
identified as a partial source in the combined provenance.

`scripts/integrate_partitioned_structure_cohort.py` checks every source receipt,
disjoint model identities, exact coverage of the original candidate FASTA, and
original/global marker-link multiplicities. It then runs the existing retained
matrix mapping and independent residue readback for the complete cohort. The
combined audit has its own partition-union status; it does not invent a single
shared prediction configuration. Six previously reused sequences remain outside
this prediction queue.

Plan: `metadata/esmfold_extended_cohort_mapping_plan.json` (389 pinned sources).
Unit: `fungal-extended-cohort-mapping-20260922.service`. Resources are one CPU
equivalent, 16 GiB memory, no swap, a 20 GiB planning output allowance and a
50 GiB free-disk gate. The 0.5–12-hour planning range begins after predecessor
validation; it is not a measured ETA. No GPU or paid resource is used.

The combined audit and model inventory are stored in
`results/predictions/extended-partition-union-audit-v1` and
`results/structures/esmfold-extended-union-v1`. Mapping and independent residue
readback will be under `results/structural_markers/esmfold-extended-complete-v1`
and `esmfold-extended-residue-readback-v1`. Controller state and logs are in
`results/structural_markers/extended-cohort-mapping-controller-v1`. The inventory union is complete; residue mapping is running and independent
residue readback follows. Native features and confidence qualification remain
subsequent requirements.

## Queued feature extraction and confidence qualification

The six-stage workflow in `scripts/advance_local_structural_features.py` is now
queued for this cohort under `fungal-extended-native-features-20260922.service`.
Its plan, `metadata/esmfold_extended_native_feature_plan.json`, binds the exact
live mapping-controller identity and requires successful full mapping and
independent residue readback before proceeding. It retains each model's source
configuration and original NPZ confidence arrays.

The stages are native Foldseek extraction, export identity/shape readback,
coordinate-feature reconstruction, lossless directional PAE export,
six-residue confidence qualification, and independent context readback against
the original arrays. The 4,363 models contain 3,924,624 residues and
3,552,697,882 PAE entries. The same previously qualified native build and source
hashes used for the other cohorts are required.

The sequential stages use at most four CPU equivalents and 16 GiB memory,
with no swap or GPU inference. The plan allows 100 GiB output and requires
150 GiB free disk before each stage. Its 0.5–36-hour runtime range is a planning
allowance after mapping finishes, not an observed ETA. Logs and state are under
`results/structural_alphabet/extended-cohort-feature-controller-v1/`; the queued
workflow is not evidence of completed feature qualification or evolutionary
analysis.
