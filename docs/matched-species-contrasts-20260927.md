# Species contrasts for the matched comparison model

All 2,786,912 selected records map to 52,675 distinct target/background node
pairs and 4,568 distinct species-weight patterns. Original record multiplicity
is retained in the pair table, while the original selected-record table remains
checksum-bound for policy, scenario and orientation. No structural qualification
or selection of favorable settings is applied at this stage.

The working species effect for a comparison is the focal species effect minus
the mean of the two background endpoint species effects:

\[
u_a - \tfrac12 u_b - \tfrac12 u_c = w_i u.
\]

This follows from an **assumed additive endpoint nuisance effect** for each
protein-pair response. It is not derived from structural-distance additivity,
and it need not be the correct evolutionary covariance. Pair interactions,
family effects, shared controls and prediction uncertainty are not captured
by this species term alone and remain necessary modeling considerations.

Weights are combined exactly when taxa coincide. For example, if the focal
species is also endpoint b, the contrast becomes one half of the focal-minus-c
effect. Endpoint order is irrelevant to these equal endpoint weights. The
original endpoint indices remain in the pair table for checking and for future
alternative models. No all-zero patterns occur in the current selected set.

The sparse matrix W has 4,568 rows and 526 columns in the fixed species-kernel
taxon order. Only 259 columns are used, and the active design has numerical rank
242. Thus 17 directions in the active species space are unidentifiable through
this design, including its constant direction. These are algebraic dimensions,
not a count of independent biological events or clades. A model must work in
the estimable space or use an explicitly specified latent covariance; the full
species design cannot simply be inverted. Structural filtering in each fitted
setting may reduce rank further.

For any of the five verified tree kernels C, a prospective species covariance
is proportional to W C W-transpose. Because each row of W sums to zero, it is
invariant to the common centering convention and to root placement in the
underlying patristic distances. For the three endpoint indices, the diagonal
quadratic form is

\[
w_i C w_i^T = \tfrac12 D_{ab}+\tfrac12 D_{ac}-\tfrac14 D_{bc}.
\]

D is patristic distance, **not squared patristic distance**. These kernel values
are in substitution-distance units, without an estimated variance multiplier.
They are not fitted variances, calibrated times or structural rates.

The producer checked the quadratic forms for all 52,675 pairs across all five
trees (263,375 checks); maximum absolute discrepancy was 6.44e-15. A separate
checker reconstructed every original selected-record multiplicity, every sparse
weight and pair assignment, and all 22,840 pattern/tree quadratic forms using
the general zero-sum distance identity. It also computed the full design rank.

Scripts: `scripts/prepare_matched_species_contrasts.py` and
`scripts/check_matched_species_contrasts.py`.
Output: `results/phylogeny/matched-species-contrasts-20260927-v1`.
Provenance: `metadata/matched_species_contrast_completed_20260927.json`.
Independent proof: `metadata/matched_species_contrast_readback_20260927.json`.

The covariance assumption, residual structure, nonlinear sequence adjustment
and uncertainty method still require validation. No phylogenetic effect has
been fitted by this preparation step, and the project remains incomplete.
