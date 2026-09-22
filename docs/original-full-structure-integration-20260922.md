# Completing integration of the original short-protein cohort

The original prediction queue and coordinate conversion contain 10,522 models.
The earlier comparative snapshot used 5,121 of these; all 5,121 exact sequence
identities occur in the completed queue. The new integration therefore adds
5,401 previously unmapped models from that queue. It does not run predictions.

The full cohort has 10,573 exact global marker/protein links spanning 95 markers
and 201 taxa. `metadata/esmfold_original_full_mapping_scope.json` records this
scope and binds the original partial mapping and complete prediction audit.
The complete mapping must replace the partial original mapping in future
unions. Combining both would double-count 5,121 models. Existing comparative
snapshots remain immutable and retain their original interpretation.

The mapping controller uses the existing
`scripts/map_audited_local_cohort.py`, followed by independent readback of
every exported retained residue against original alignment and prediction
evidence. Its frozen plan is `metadata/esmfold_original_full_mapping_plan.json`;
the service is `fungal-original-full-mapping-20260922.service`.
Outputs are under `results/structural_markers/esmfold-original-complete-v1`
and `results/structural_markers/esmfold-original-residue-readback-v1`.

The pre-launch resource allowance is one CPU equivalent, 16 GiB RAM, no swap,
20 GiB output, a 50 GiB free-disk gate, and 0.5–12 hours of uncalibrated runtime
planning. Mapping has started; completion requires its bound final receipt.

Native structural-feature extraction and full confidence qualification are
queued behind the exact mapping process via the existing
`scripts/advance_local_structural_features.py`. This repeats the established
six stages: native extraction, export identity readback, independent coordinate
feature reconstruction, PAE export, confidence qualification, and original
PAE-context readback. The plan is
`metadata/esmfold_original_full_native_feature_plan.json`, with service
`fungal-original-full-native-features-20260922.service`.

The feature stage covers 3,355,507 protein residues and 1,202,825,245 PAE matrix
entries. Its allowance is four CPU equivalents, 16 GiB RAM, no swap, 100 GiB
output, a 150 GiB free-disk gate and 0.5–24 hours of uncalibrated planning after
mapping completes. Original coordinates and prediction configuration remain
bound to each model. No GPU prediction or paid resource is involved.

Subsequent work must integrate the complete original cohort with the follow-on,
ecology, long and extended cohorts, retain configuration differences, regenerate
paired inference inputs, and rerun the relevant phylogenetic analyses. A larger
model count alone does not establish increased qualified alignment coverage or
support for any evolutionary hypothesis.
