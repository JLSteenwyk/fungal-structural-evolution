# Reference alignment numerical diagnostics

The whole-protein reference batch completed 130,164 directed dispositions:
125,836 aligned and 4,328 unavailable under the pLDDT70 input mask. Its original
strict numeric audit failed at a reverse pLDDT70 alignment of two residues.
Recomputed RMSD 0.215449024570517 Å differs from printed 0.21 Å by
0.005449024570517 Å, beyond the unchanged 0.00501 Å tolerance. The two centered
coordinate sets each have rank one and do not uniquely specify a rigid rotation.
This local finding does not establish why the native implementation printed its
value or show that this is the only discrepancy.

The original failed audit remains intact. Separate stages now cover the full
reference batch:

1. `scripts/diagnose_reference_alignment_rmsds.py` checks every disposition,
   source mapping and successful numeric result, recording all RMSD failures
   rather than aborting at the first. The result explicitly lacks scientific
   acceptance and retains discrepancy flags.
2. `scripts/assess_reference_alignment_geometry.py` waits for the identified
   diagnostic producer and a complete receipt. It rechecks every successful
   residue mapping and numeric record, then calculates coordinate rank, singular
   values and curvature of the optimal proper rotation.
3. `scripts/readback_reference_alignment_geometry.py` waits for the complete
   geometry result, independently rebuilds residue indexes, checks singular
   values with LAPACK gesvd and tests curvature using a quaternion eigensystem.

The geometric calculations preserve the previously tested domain routines;
separate variants handle the reference cohort's two model/input bundles and
model/version/mask identities. They do not change the running or historical
scripts. Each stage has one CPU, 16 GiB memory and no swap. Geometry and readback
allow 0.2–8 hours each and 2/1 GiB output respectively, without new predictions.
Both waiting processes have been verified live; full geometry results and their
verification remain pending. Numerical uniqueness is not prediction accuracy.

Plans and launch identities are in
`metadata/duplication_reference_{rmsd_diagnostic,geometry,geometry_readback}_{plan,launch}_20260927.json`.
The first failure is archived in
`metadata/duplication_reference_first_rmsd_failure_20260927.json`.
Downstream scientific eligibility and handling of short/degenerate mappings
remain unresolved. The original strict order-summary gate remains unsatisfied.
