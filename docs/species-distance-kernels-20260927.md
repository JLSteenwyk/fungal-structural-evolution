# Audited species-distance kernels for model sensitivity

Prepared full 526-tip patristic-distance matrices for five completed, audited
tree alternatives: profile and MAFFT guide trees and the profile/profile,
profile/MAFFT and MAFFT/profile PMSF runs. The fourth PMSF run and hybrid-exclusion
sensitivities remain unfinished; these five trees are not declared the final
species framework. All branch lengths remain in inferred substitution units.

For each tree, let X contain root-to-tip edge indicators multiplied by square
roots of edge lengths. Then D_ij = ||X_i - X_j||² is the patristic distance.
With H = I - 11'/n, the centered positive-semidefinite kernel is
C = (HX)(HX)' = -HDH/2. D is not squared again. The producer computes both forms
and checks eigenvalues at numerical tolerance. Kernel centering is over the
complete shared 526-tip grid; it is not an estimate of trait covariance.

Independent verification reconstructed every one of the 690,375 unordered
pair/tree distances using DendroPy instead of the producer's Bio.Phylo path
features. It then explicitly rerooted each tree on a terminal edge and checked
all distances again. All kernel entries, recovered distances and eigenvalues
also passed alternative linear-algebra checks. Maximum distance disagreement
was 4.441e-15 and maximum reroot disagreement 2.665e-15. Near-zero negative
eigenvalues (about 1e-15) are numerical roundoff, not evidence of negative
biological variance; no inverse kernel or automatic regularization is exported.

For a future zero-sum endpoint-contrast design W, W C W' supplies a possible
root-independent dependence kernel. Choosing W and fitting its variance scale
would still be explicit modeling assumptions; structural distances are not
Brownian traits by virtue of this construction. This is not time calibration,
validation of an evolutionary covariance process, or a correction already
applied to the matched comparisons. Rooting, model adequacy and tree sensitivity
remain separate scientific questions.

All 362 taxa in the full target/background matching graph occur in the common
526-tip grid. Every graph node was checked, including currently unselected and
structurally unavailable records. Proof:
`metadata/species_kernel_graph_taxon_coverage_20260927.json`.

Plan: `metadata/species_distance_kernel_plan_20260927.json`.
Output: `results/phylogeny/species-distance-kernels-20260927-v1`.
Full proof: `metadata/species_distance_kernel_readback_20260927.json`.
Scripts: `prepare_species_distance_kernels.py`,
`readback_species_distance_kernels.py`, and `check_species_kernel_graph_taxa.py`.
NPZ artifacts preserve all distances, centered kernels and eigenvalues; the
common ordered labels are in `taxa.json`. Source trees/audits and every output
are checksum-bound. The producer ran under one CPU, 8 GiB and no swap; the
prelaunch allowance was 1 GiB output and 1–30 minutes. Independent verification
used one BLAS thread. No new tree inference, GPU prediction or paid resource.
