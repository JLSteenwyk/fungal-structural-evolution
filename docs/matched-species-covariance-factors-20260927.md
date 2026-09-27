# Compact factors for the working species covariance

The complete species contrast covariance now has verified factors for all five
tree alternatives. Each factor has 4,568 pattern rows and 242 columns. These
provide the prospective phylogenetic term without storing a dense covariance
for millions of selected records or inverting an unidentifiable species design.

For the full contrast design W, singular value decomposition gives W = Q R
on its numerical rank-242 space. Q has orthonormal columns and R maps the 526
species columns to that space. For each verified species kernel C, the reduced
matrix B = R C R-transpose has a successful Cholesky decomposition B = L L-transpose.
The exported factor F = Q L therefore satisfies F F-transpose = W C W-transpose.
All five reduced matrices are positive definite; the smallest eigenvalues range
from approximately 0.000486 to 0.000642. The computation required no diagonal
jitter, covariance eigenvalue clipping or inverse of the species matrix.

The numerical rank threshold is max(W dimensions) times machine epsilon times
the largest singular value (approximately 1.19e-11). Reconstructing W from the
retained factors was checked against every original design entry. Dropped
singular directions are numerical null directions, not selected biological
effects. Each tree's full pattern covariance was then reconstructed and checked
against W C W-transpose in memory-bounded row blocks.

An independent verifier parses the original trees with DendroPy and constructs
tip-to-edge incidence features weighted by square root of branch length. Since
the contrast weights sum to zero, contrasts of these features yield the same
covariance without centering and without reading C for the numeric comparison.
It checked all 104,333,120 ordered pattern-pair/tree entries, with maximum
absolute discrepancy 2.49e-14, as well as the complete basis, species projection,
reduced matrices and Cholesky factors. This establishes numerical equivalence
for the specified model input, not biological adequacy of that covariance.

For a fitted setting, select factor rows using each retained record's
`species_pattern_id`. Repeated patterns must retain repeated observations and
their control/family identifiers. A fitted variance multiplier and positive
residual component are still needed. Background and family dependence,
nonlinear sequence adjustment, uncertainty and model diagnostics remain
separate requirements. Qualification can further reduce the rank of the
selected rows; the exported full-design factor does not assert otherwise.

Scripts: `scripts/factor_matched_species_covariance.py` and
`scripts/readback_matched_species_factors.py`.
Output: `results/phylogeny/matched-species-covariance-factors-20260927-v1`.
Provenance: `metadata/matched_species_factor_completed_20260927.json`.
Independent proof: `metadata/matched_species_factor_readback_20260927.json`.
Both commands completed successfully with numerical-library threads restricted
to one. No GPU inference or paid resources were used. No evolutionary effect
has been fitted by this factorization.
