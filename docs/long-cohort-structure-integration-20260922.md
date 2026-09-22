# Integrating the completed 513–768-residue marker cohort

**September 22 update:** mapping and independent readback passed for all 5,510
models, 5,544 marker links and 1,984,383 retained residue links. Native feature
extraction and export readback also passed; coordinate-feature reconstruction
is running, with PAE qualification and context readback still pending. Receipts:
`metadata/esmfold_long_cohort_mapping_completed.json` and
`metadata/esmfold_long_cohort_residue_readback.json`. The launch descriptions
below retain their original planning context.

All 5,510 models in this length cohort passed the independent prediction
artifact audit and sequence-explicit mmCIF conversion. They match 5,544 global
marker-protein records across 108 markers and 279 taxa by exact complete
amino-acid sequence. The originating prediction links and the full global
exact-sequence link set agree for this cohort; no matching taxon record is
discarded through sequence deduplication.

On September 22, `fungal-long-cohort-mapping-20260922.service` began mapping
these models to retained columns of the frozen profile alignment matrix. The
controller `scripts/map_audited_local_cohort.py` checks the prediction audit,
conversion receipt, inventory identity set, and global link multiplicities.
It then runs the existing `map_marker_structures.py` followed by
`audit_local_residue_mapping_global_links.py`.

The independent readback reconstructs every expected non-gap retained residue
from the source Stockholm alignments, checks the corresponding matrix amino
acid and model position, and compares confidence with the original prediction
NPZ arrays at the PDB serialization precision. Model, prediction and coordinate
hashes remain bound to their original provenance. This validates mapping; it
does not independently infer the alignment or establish model accuracy.

The plan `metadata/esmfold_long_cohort_mapping_plan.json` pins 392 script and
data sources, including all relevant unaligned marker FASTAs, Stockholm
alignments, alignment receipts, the matrix, and the audited local inventory.
The service uses one CPU equivalent, a 16 GiB hard memory cap (8 GiB planning
allowance), no swap, a 20 GiB output allowance and a 50 GiB free-disk gate.
The planning runtime allowance is 0.5–8 hours. No GPU or paid resource is used.

Mapping outputs are under
`results/structural_markers/esmfold-long-complete-v1`, and the independent
readback is under `results/structural_markers/esmfold-long-residue-readback-v1`.
Controller logs and the exact global link list are under
`results/structural_markers/long-cohort-mapping-controller-v1`.
At launch, mapping is running and independent readback is queued immediately
after it; completion is unproven.

PAE binding, native structural-alphabet extraction, coordinate-feature checks,
confidence qualification, integration with the newer 769–1,024-residue models,
and evolutionary analyses remain subsequent steps. The previously published
89-marker ESMFold results do not yet incorporate this cohort.

The six-stage native-feature follow-up is now queued under
`fungal-long-native-features-20260922.service`, using
`scripts/advance_local_structural_features.py` and
`metadata/esmfold_long_native_feature_plan.json`. It waits for the exact live
mapping-controller identity, then requires the completed controller receipt
and independent residue-readback receipt before any feature work begins.

The stages use the existing qualified implementations: pinned Foldseek 3Di
extraction; native export identity/shape readback; coordinate reconstruction
of spatial partners and ten descriptors; lossless directional PAE export;
six-residue confidence qualification; and independent context readback against
the original local NPZ arrays. The last check evaluates all 36 ordered residue
pairs in each valid six-residue context. It also verifies that coordinate
arrays and their summary counts remain unchanged. Each stage must retain the
same 5,510-model universe and mapping provenance.

This cohort contains 3,386,360 residues and 2,109,723,548 PAE entries. The
sequential follow-up has a four-CPU-equivalent cap, 16 GiB memory, no swap,
100 GiB planning output allowance, and a 150 GiB free-disk gate before each
stage. The 0.5–24-hour planning range is not a measured ETA. No GPU prediction
or paid resource is used. Native encoder source files are retrieved from the
same immutable public commit and compared with the previously qualified source
hashes; a different build cannot silently enter this cohort.

Controller state, stage logs and eventual completion evidence are under
`results/structural_alphabet/long-cohort-feature-controller-v1/`. Queuing these
stages does not establish completion or validate biological structural changes.
