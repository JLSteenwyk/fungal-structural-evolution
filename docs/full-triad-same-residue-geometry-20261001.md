# Full same-residue duplicate/reference geometry — October 1, 2026

The full geometric comparison and independent quaternion reconstruction have
completed on existing CPUs. Closure verified 464,600 source/artifact bindings
and both exact original producer/reader journals. The production source is
the independently verified full mapping grid
described in [the correspondence workflow](full-triad-common-residues-20260930.md).

All 27,056 source-ready ordered physical triples are included under two masks,
eight native AB/AR/BR order combinations and two correspondence definitions:
**432,896 mapping states and 865,792 fit dispositions**. No favorable order or
coverage outcome is used to select the scheduled cohort. The immutable upstream
design retains all 283,409 contexts, 214,461 reference ties and 428,922 logical
sides, including unscheduled and missing states.

Full-data preflight exhausted the entire grid and rechecked its source bindings.
It found 713,512 dispositions with eligible coordinate inputs, 105,732 short-core
dispositions and 46,548 source-excluded dispositions. These were prefit statuses:
all 713,512 eligible dispositions subsequently had unique fits at the declared
numeric tolerance. The short and source-excluded dispositions remain explicit.
Within a physical
triple and mask, order/core alternatives repeat some residue sets. There are
77,805 distinct eligible sets, totaling 16,923,795 residue occurrences and
233,415 proper pair fits per implementation after identical-core reuse. The
142,847 needed PDB inputs total 3,841,629,229 bytes. Preflight closure checked
321,735 bindings and its exact original process journal. This is an exhaustive
full-data work estimate, not a pilot or a measured geometry runtime.

For each correspondence definition, all three rigid fits use the identical
original-position residue triples. The producer reuses the pinned proper SVD
fit and geometry routines; the reader independently solves the quaternion
eigenproblem. Every distance, signed AR-minus-BR contrast, sequence identity,
confidence summary, rotation classification and output row is reconstructed.
The signed contrast retains original A/B orientation and is a descriptive
distance difference, not a directional evolutionary rate. The reader also
checks the common-core RMSD metric bound.

Written PDB checksums, residue numbers, amino-acid identities, coordinate
finiteness and confidence values are checked before fitting. Both original and
retained input lengths are exported, with retained coverage explicitly labeled.
Core coverage uses original full-protein lengths for all three roles under
both masks. The six fixed screens remain 30/50 residues crossed with
50%/70%/90% original coverage. Shared-core pass flags, inherited both-order
three-edge pass flags and their conjunction are exported separately. Inherited
exclusions retain the responsible edge. Failed sources and cores with fewer
than three residues have blank geometric metrics; degenerate fits retain their
computed distances and remain excluded.

Synthetic fixtures passed 48 mapping states and 96 fit rows, including rigid
transforms, reflection, collinearity, unequal protein lengths, nonconsecutive
confidence masks, short cores and inherited source/coverage exclusions. The
independent maximum distance/contrast difference was approximately 3.02e-14.
All 13 rehashed false exports were rejected, including wrong denominators,
changed residue hashes or model roles, promoted degenerate/short cores, cleared
source exclusions and favorable inherited screens. Proof I/O was stubbed only
in these software fixtures; full production uses the closed actual-data gate.

The full reader reproduced every row, distance, signed contrast, confidence,
identity and screening decision. Maximum absolute RMSD/contrast disagreement
was 5.107e-14 Å. At 50 residues/70% original coverage, combined core/three-edge
passes number 66,984 reference-common versus 59,480 cycle-consistent full-mask
states, and 45,440 versus 43,624 pLDDT70 states. These counts include eight-order
alternatives and are not independent events, physical triples or calibrated
effect estimates. Full logical-context linkage and order/mask robustness
qualification remain ahead.

Each serial stage had two CPU cores, 32 GiB RAM, no swap and one BLAS thread.
The output allowance is 32 GiB with a 100 GiB disk reserve. The planning range
of 0.5–12 hours per stage remains uncalibrated. No GPU or paid resources are used.

Reproducibility records:

- [Closed full-data work preflight](../metadata/full_triad_same_residue_fit_preflight_verified_20261001.json)
- [Measured work and resource estimate](../metadata/full_triad_same_residue_fit_actual_resources_20261001.json)
- [Software fixture results](../metadata/full_triad_same_residue_fit_fixture_validation_20261001_v2.json)
- [Immutable production plan](../metadata/full_triad_same_residue_fit_plan_20261001.json)
- [Original producer, reader and closure handles](../metadata/full_triad_same_residue_fit_pipeline_started_20261001.json)
- [Verified full geometry completion](../metadata/full_triad_same_residue_fit_completed_20261001.json)
- [All 58 unit-labeled fit/screen count rows](tables/full_triad_geometry_dispositions_20261001.tsv)

Large tables, source hashes and independent readback remain outside Git in
`results/structural_comparisons/full-triad-same-residue-fits-20261001-v1/`.
The small completion record is
`metadata/full_triad_same_residue_fit_completed_20261001.json`, with status
`complete_verified_full_triad_same_residue_geometry`.

Logical-context linkage, native assignment and reference-choice sensitivities,
sequence-locked/domain/PAE/orientation controls, prediction calibration,
phylogenetic dependence and calibrated inference remain required. The fourth
crossed PMSF run separately passed its full profile/tree/bootstrap audit and
[37-binding, two-journal closure](../metadata/pmsf_fourth_readback_completed_20261001.json).
All-four topology/consensus comparison is now independently closed in
[the species-tree sensitivity workflow](pmsf-four-run-sensitivity-20261001.md).
Rooting, model adequacy and taxon/marker/hybrid sensitivities remain required
before accepting a species framework. All eight
scientific aims remain incomplete; GPU predictions remain paused.
