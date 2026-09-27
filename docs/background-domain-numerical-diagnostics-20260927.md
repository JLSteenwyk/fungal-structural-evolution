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

A separate full diagnostic is now running:
`scripts/diagnose_background_domain_rmsds.py --plan
metadata/background_domain_rmsd_diagnostic_plan_20260927.json`.
It preserves all source, input, mapping and identity checks, but records every
RMSD discrepancy explicitly instead of aborting at the first. The original
failed strict audit and its tolerance remain unchanged. The new diagnostic
cannot establish scientific acceptance. Output:
`results/structural_comparisons/background-domain-rmsd-diagnostic-20260927-v1`.

The diagnostic was confirmed live under one CPU, 16 GiB and no swap; budget is
0.2–8 hours and 2 GiB output, using existing resources without predictions or
new charges. Full diagnostic readback, geometry/rotation qualification and
matched-control inference remain downstream. Both target and background
comparisons must undergo compatible numerical and coverage criteria.
