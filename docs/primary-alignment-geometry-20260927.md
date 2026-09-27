# Primary whole-protein alignment geometry

The main duplicate-pair comparisons require a check of whether each aligned
coordinate set determines a unique rigid rotation. A small RMSD alone does not
establish this: very short or collinear mappings can leave rotations ambiguous.
The new full-cohort stage extends the completed reference-geometry method to
all successful primary alignments, preserving both input orders and confidence
masks. No outcome-based selection of pairs is used.

The controller waits for the recorded strict primary numeric auditor's exact
process identity and successful terminal service state. It requires the audit
receipt to bind the completed producer and its native checkpoint manifest.
A failed strict audit stops this stage; the controller does not bypass or loosen
RMSD checks. The older failed reference audit and its separately documented
short-alignment diagnostics remain unchanged.

For each successful primary checkpoint, the stage rechecks native-text mapping,
sequence identity and RMSD against hashed materialized PDBs and the strict
numeric table. It computes coordinate rank, singular spectra and proper-rotation
curvature, then checks these with the existing alternate-SVD and quaternion
implementation. Every successful checkpoint must appear exactly once. Counts
of excluded native dispositions remain tied to the upstream complete audit.
Near-zero quaternion gaps are explicitly counted; cross-method agreement does
not prove rank at machine precision or stability to coordinate uncertainty.

Output: `results/structural_comparisons/primary-alignment-geometry-20260927-v1`.
Run `scripts/assess_primary_alignment_geometry.py --plan
metadata/primary_alignment_geometry_plan_20260927.json` from the repository root.
The pinned plan and launch record preserve inputs, scripts and dependency
identity. Resources: one CPU with BLAS/OpenMP restricted to one thread,
16 GiB memory, no swap, and 2 GiB planned output. The 0.5–12 hour planning
allowance excludes the primary producer and audit wait; it is not a measured ETA.
No GPUs or paid infrastructure are used.

## Integration checks and current status

The native mapping and geometry path passed on 32 saved successful primary
checkpoints. Each geometry row was serialized to TSV and read back through the
quaternion checker. Four collinear/degenerate coordinate fixtures of lengths
1, 2, 3 and 12 retained the expected degeneracy classifications. Exact checked
native paths/hashes are in
`metadata/primary_alignment_geometry_integration_checks_20260927.json`.
These integration checks do not establish full-cohort completion.

The controller is queued behind the running primary audit; its identity is in
`metadata/primary_alignment_geometry_launch_20260927.json`. Full geometry output,
a complete serialized-table readback and the downstream shared-residue triad
comparisons remain pending. This stage tests numerical identifiability, not
prediction accuracy, structural asymmetry or ancestral change.
