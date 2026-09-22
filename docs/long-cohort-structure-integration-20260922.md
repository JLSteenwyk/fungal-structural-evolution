# Integrating the completed 513–768-residue marker cohort

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
