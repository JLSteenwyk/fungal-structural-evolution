# Background-domain numerical diagnostics

The full 66,929-domain-pair batch completed successfully: 267,716 directed
mask/order outcomes, comprising 267,246 aligned and 470 input-unavailable cases.
Its strict independent numeric audit then failed at a pLDDT70 mapping of two
residues. The exact checkpoint and source hashes are retained in
`metadata/background_domain_first_rmsd_failure_20260927.json`.

The first failed mapping has recomputed RMSD 0.147182554481412 Å versus native
printed 0.14 Å, discrepancy 0.007182554481412 Å > unchanged 0.00501 Å tolerance.
Both centered coordinate sets have rank one and the rotation is nonunique.
The analytic two-point formula (half the absolute difference between interpoint
distances) independently gives 0.147182554481411 Å. This confirms the coordinate
minimum but does not explain the native implementation's reported value or
prove that other rows are unaffected.

The separate full diagnostic completed successfully:
`scripts/diagnose_background_domain_rmsds.py --plan
metadata/background_domain_rmsd_diagnostic_plan_20260927.json`.
It preserves all source, input, mapping and identity checks, but records every
RMSD discrepancy explicitly instead of aborting at the first. The original
failed strict audit and its tolerance remain unchanged. The new diagnostic
cannot establish scientific acceptance. Output:
`results/structural_comparisons/background-domain-rmsd-diagnostic-20260927-v1`.

The completed diagnostic reconstructs all 267,246 successful mappings. Five
RMSDs exceed the original tolerance, all involving two aligned residues. In
all, 13 alignments contain two residues and 267,233 contain at least three.
The completion check verifies every successful checkpoint key, all artifact
hashes, counts, classifications, and all pinned sources; it is recorded in
`metadata/background_domain_rmsd_diagnostic_completed_20260927.json`.
The original strict audit remains failed, and scientific eligibility remains
unestablished.

Full geometry assessment is running under
`metadata/background_domain_geometry_plan_20260927.json`. An independent
quaternion/alternative-SVD readback is queued under
`metadata/background_domain_geometry_readback_plan_20260927.json` and waits
for the exact producer process identity. Each stage uses one CPU, 16 GiB,
no swap, and a planning estimate of 0.2–8 hours; output budgets are 2 GiB and
1 GiB, respectively. Mathematical functions match the existing domain
implementation exactly. Point, two-point, planar, reflected and random
fixtures pass; corrupted rank, curvature, status and ratio are rejected.
These stages assess numerical identifiability, not prediction uncertainty.

Matched-control inference remains downstream. Both target and background
comparisons must undergo compatible numerical and coverage criteria.


All 13 two-residue mappings were independently checked using half the absolute
difference of their two inter-residue distances. This verifies the minimum
RMSD without matrix fitting; maximum disagreement with the diagnostic is
1.027e-15 Å. All 13 lack a unique rotation, including all five RMSD discrepancies.
Full results and hashes are retained in
`metadata/background_domain_short_geometry_readback_20260927.json`.

The numerical order summary is queued under
`metadata/background_domain_usable_orders_plan_20260927.json`. After full
geometry verification, it will retain all 66,929 pairs and both masks/orders,
blanking metrics for RMSD discrepancies, fewer than three paired residues,
or numerically nonunique rotations. An independent full table readback runs
automatically afterward. The exclusion and summary functions are unchanged
from the previously verified reference implementation. This preserves
comparability without declaring coverage or prediction-confidence eligibility.
The summary and verification use one CPU, 16 GiB, no swap, up to 2 GiB output,
and an estimated 0.1–2 hours after dependencies finish. No GPU predictions.


## Completed geometry and numerical order verification

The independent geometry verifier passed all 267,246 directed alignments.
Exactly 13 fits are degenerate, identical to the complete analytic two-residue
set; maximum scaled quaternion-curvature disagreement is 2.525e-15. Both
producer and verifier terminated successfully. Completion is archived in
`metadata/background_domain_geometry_completed_20260927.json`.

The order summary also passed full independent verification: all 133,858
pair/mask rows and 3,206,784 numeric values were checked. There are 267,233
numerically usable directions and 13 excluded directions; 66,929 full-mask
pairs and 66,686 pLDDT70 pairs have both orders usable. Three pLDDT70 pairs
have only one usable direction, and 240 have neither. Proof:
`metadata/background_domain_usable_orders_readback_20260927.json`.
Coverage and geometry qualification for both target and background cohorts
has also completed; see `docs/domain-coverage-geometry-20260927.md`.
These numerical checks do not establish prediction accuracy or a duplication effect.
