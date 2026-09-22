# Completed extended marker predictions and validation

The 769–1,024-residue prediction queue finished on September 22, 2026.
Its 4,363 candidate sequences are present in three disjoint output directories:
2,318 in the original interrupted run, 1,022 in the GPU 0 continuation, and
1,023 in the GPU 1 continuation. The original interruption receipt remains
unchanged; its remaining 2,045 sequences were completed by the two continuations.
The separate 513–768-residue cohort contains 5,510 audited, converted predictions.
These are marker cohorts, not predictions for the full fungal proteomes.

Independent validation of the 4,363 new predictions is **running**, not yet
complete. The CPU-only controller is
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
