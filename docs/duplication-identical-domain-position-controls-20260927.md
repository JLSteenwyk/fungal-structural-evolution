# Identical-domain sequence-position controls — September 27, 2026

Completed direct coordinate comparisons for all 48 robust candidates whose
complete domain interval strings are identical across every alternative boundary.
The analysis includes all eight candidates previously identified as sensitive
to residue correspondence; it does not select only those eight.

The 48 candidates resolve to 75 distinct oriented interval pairs. Each was fit
using exact sequence offsets under two masks: all interval CA positions, and
the intersection of positions with pLDDT at least 70 in both predictions.
No structural alignment search chooses the paired residues. All 150 fits have
unique rotations at the numerical tolerance and pass the descriptive minimum
of 30 paired residues and 70% original-interval coverage.

| Mask | Interval pairs | Minimum RMSD (Å) | Median RMSD (Å) | Maximum RMSD (Å) |
| --- | ---: | ---: | ---: | ---: |
| All positions | 75 | 0.153 | 0.512 | 4.211 |
| Both pLDDT ≥70 | 75 | 0.153 | 0.375 | 3.421 |

The 75 interval pairs are boundary alternatives within 48 candidate combinations,
not independent evolutionary replicates. Masking changes the fitted residues;
the reduction in median distance is descriptive and does not calibrate an error
bound. pLDDT is retained as a model-quality measure, not evidence that any predicted
shape is experimentally correct.

These fits show predicted coordinate differences after exact sequence-position
pairing. They do **not** validate the original three-protein directional contrast.
The next correspondence control must carry identical duplicate positions into
both duplicate/reference comparisons on a shared reference core, retaining
reference and confidence alternatives. Full-chain sequence/context differences
and prediction uncertainty remain possible explanations for the pair distances.
No ancestral change, selection or biological sequence-independent divergence
is established.

## Reproduction and verification

Run `scripts/fit_identical_domain_positions.py --plan
metadata/duplication_identical_domain_position_plan_20260927.json`, followed by
`scripts/check_identical_domain_positions.py --plan
metadata/duplication_identical_domain_position_plan_20260927.json --output
<new-readback-path.json>`.

The immutable output directory
`results/structural_comparisons/duplication-identical-domain-position-fits-20260927-v1`
contains `position_fits.tsv`, `candidate_pair_links.tsv` and a hash-bound receipt.
Every sequence offset, full interval identity, mask, mean confidence and numerical
geometry status is retained. Candidate links preserve role and relationship
class, including the earlier correspondence flag.

The independent checker read all 150 PDBs, reconstructed sequences and offsets,
and used quaternion rotations to verify the producer's SVD fits. It checked the
complete candidate/pair/mask grid, all RMSDs, geometry, confidence and screen
values. Maximum RMSD disagreement was 2.18 × 10⁻¹⁴ Å. See
[completion evidence](../metadata/duplication_identical_domain_position_completed_20260927.json).
The resource plan allowed one CPU, 2 GiB memory and 1 GiB output. No predictions
or native alignment searches ran; the GPU pause remains in force.
