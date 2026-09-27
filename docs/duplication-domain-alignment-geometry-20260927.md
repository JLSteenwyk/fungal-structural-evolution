# Domain alignment geometry assessment

The geometry assessment covers all 280,824 successful directed domain
alignments, not only pairs passing a coverage screen. Unavailable-input
comparisons remain in the upstream disposition ledger. The separate numerical
diagnostic, original checkpoint hashes and PDB manifests bind this analysis to
the same data. Residue mappings and all diagnostic metrics are checked again
before geometry is calculated; the existing native RMSD discrepancy remains
explicit rather than being corrected or accepted.

For centered matched coordinate matrices A and B, the output records each
coordinate singular spectrum as RMS widths in ångströms (singular values divided
by the square root of the matched count). Numerical ranks use float64 epsilon
times max(n,3) times the largest coordinate singular value. Ratios of second
to first widths are reported continuously to expose near-collinearity; no
post hoc cutoff is used to label predictions reliable.

Coordinate rank alone does not establish a unique optimal proper rotation.
For AᵀB = U diag(s1,s2,s3) Vᵀ with ordered nonnegative singular values, define
δ = sign(det(UVᵀ)). The proper-rotation solution uses U diag(1,1,δ) Vᵀ.
The squared-error rotation Hessian has eigenvalues twice the pairwise sums of
(s1,s2,δs3). Thus its smallest eigenvalue is 2(s2+δs3). The output records
`s2+δs3`, its ratio to s1, the determinant correction and all cross singular
values. A fit is `unique_at_numeric_tolerance` only when this minimum curvature
exceeds eps(float64) × max(n,3) × s1. Otherwise it is numerically degenerate.
This identifies a local flat rotational direction; it does not prove robustness
to prediction error or biologically relevant perturbations.

The tests cover one-point and collinear degeneracy, a noncollinear planar
case, isotropic reflected correspondences, swapped input orders, and independent
rigid coordinate transformations. For a nontrivial reflected/noisy example,
an independently evaluated finite-difference Hessian agrees with the curvature
formula. Planarity alone is not rejected; a noncollinear planar correspondence
can determine a unique proper rotation.

The [pinned plan](../metadata/duplication_domain_alignment_geometry_plan_20260927.json)
allows one CPU, 16 GiB memory, no swap, and 1 GiB output. The prior all-record
readback took about five minutes; 5–30 minutes is a planning range for this
expanded assessment, not a guaranteed completion time. This task runs no
structure prediction or native alignment and incurs no external charge.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/check_domain_alignment_geometry_cases.py
OPENBLAS_NUM_THREADS=1 python scripts/assess_domain_alignment_geometry.py --plan metadata/duplication_domain_alignment_geometry_plan_20260927.json
```

Existing outputs are never overwritten. A completed producer receipt still
requires independent readback before these measurements are incorporated into
the domain-triad dataset. Numerical identifiability is a technical property,
not a substitute for confidence, domain-boundary, sequence-evolution or
phylogenetic controls.

## Independent full-cohort readback

The [readback plan](../metadata/duplication_domain_geometry_readback_plan_20260927.json)
waits for the exact producer PID/creation time/command to finish and requires a
completed, hash-bound geometry receipt. It checks the exact successful-alignment
key set and reconstructs every matched coordinate pair from hashed PDBs and
native alignment strings. Original scripts and artifacts remain unchanged.

Coordinate and cross spectra are recomputed using LAPACK `gesvd`, rather than
the producer's NumPy SVD implementation. A separate 4×4 symmetric quaternion
matrix supplies the optimal rotation objective and curvature: half the gap
between its largest two eigenvalues equals the producer's minimum curvature.
This uses a different formulation from the determinant-corrected 3×3 SVD.
The scaled cross-method agreement tolerance is 1e-10; it does not certify
machine-precision rank. Near-zero quaternion gaps are counted separately.
Numerical rank, recorded algebra, ratios and original tolerance classifications
are also checked explicitly. No tolerance in the native RMSD diagnostic changes.

Tests cover point, two-point, planar, reflected and random coordinates. Altered
curvature, rank, status and relative-curvature fields must fail the checker.
The readback uses one CPU and at most 16 GiB RAM with no swap; it does not rerun
USalign or use GPUs. Its output remains pending until the full scan completes.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/check_domain_geometry_readback_cases.py
OPENBLAS_NUM_THREADS=1 python scripts/readback_domain_alignment_geometry.py --plan metadata/duplication_domain_geometry_readback_plan_20260927.json
```
