# Domain-anchored coordinate displacement

For the 13 opposing whole-protein/domain cases, the coordinate analysis uses
whole-protein residue correspondences and partitions each common triple by the
three annotated domain intervals. Triples whose three residues are inside form
the anchor; triples whose three residues are outside form the evaluation set.
Mixed triples are exported separately and excluded from these two fits. This
preserves one correspondence definition throughout each measurement.

The analysis retains 26 interval triads, both confidence masks, eight alignment
order combinations and two mapping definitions: 832 partitions. Each partition
is assessed for A–B, A–reference and B–reference, producing 2,496 pair records.
Domain anchors require at least three residues and unique numerical rotation
geometry. Anchor counts and curvature are retained; this numerical requirement
is not a scientific confidence or length screen.

For each usable anchor, a proper rigid transformation is fitted to the inside
residues and applied unchanged to outside residues. The outside residual is also
compared with an independent fit to the same outside coordinates when at least
three outside residues exist. This latter optimum cannot exceed the residual
under the domain transform; the inequality is checked. Small or missing outside
sets and degenerate anchors are retained explicitly.

A separate quaternion fit reproduces both anchor and outside residuals. Synthetic
fixtures test a pure rigid transform and translation of the outside region while
the domain core stays fixed. Partition membership is checked independently with
Boolean arrays; complete residue partitions and pair measurements are exported
and read back.

These measurements do not identify rigid interdomain rotation by themselves.
Outside residues can include multiple domains, flexible linkers, termini and
regions with uncertain correspondence. Domain orientation claims require a
second confidently matched domain, appropriate uncertainty checks and independent
prediction evidence. The analysis also uses whole-protein mappings rather than
native domain-alignment mappings; disagreement between those mappings remains a
separate sensitivity question.

Run `scripts/measure_domain_anchored_displacement.py` with one CPU and one BLAS
thread. Outputs are in
`results/structural_comparisons/domain-anchored-displacement-20260927-v1`.
The [case dataset](whole-domain-inspection-cases-20260927.md) supplies all chosen
references and cases. No new structure prediction is required.

The run completed all 832 partitions and 2,496 pair records with usable numerical
anchor geometry and outside residues. Independent quaternion residuals agreed
within 1.36e-13 Å; synthetic fixtures and full serialization checks passed.
Anchor sizes range from 24 to 320 residues; outside sets contain 82–781 residues,
and mixed sets contain 0–16 triples. Thus, numerical fit completion must not be
reported as passing the original domain coverage screens: whole-protein mappings
can yield different anchor sets, including anchors below 30 residues. Original
interval coverage and anchor-size screens must be reapplied to these new sets
before selecting scientifically interpretable comparisons.

The completion receipt is
`metadata/domain_anchored_displacement_completed_20260927.json`. Biological
localization, confidence qualification and prediction-error calibration remain
separate stages.

## Original-length coverage checks

`screen_domain_anchored_displacement.py` completed all 4,992 partition/screen
rows and 78 case/screen rows. The anchor denominator is each original annotated
domain interval; the outside denominator is the original protein length minus
that interval. Every one of the three proteins must meet the threshold. Outside
screening is a separate sensitivity requirement and does not change the anchor
flag. Missing or confidence-masked residues remain in the original denominators.
Integer comparisons and independent rational arithmetic agreed on every decision;
partition coordinates and all serialized outputs were checked.

| Minimum residues / coverage | All anchors pass, of 13 cases | Both anchor and outside pass, of 13 cases |
| --- | ---: | ---: |
| 30 / 50% | 11 | 8 |
| 30 / 70% | 11 | 3 |
| 30 / 90% | 6 | 1 |
| 50 / 50% | 11 | 8 |
| 50 / 70% | 11 | 3 |
| 50 / 90% | 6 | 1 |

The all-alternative flags require every interval definition, reference, mask,
order and mapping definition, not merely one passing fit. At n50/c70, the three
cases passing both requirements are Heliocybe OG0000054/PF00043.31, Jaapia
OG0000294/PF00009.34 and Phycomyces OG0002812/PF08267.19. These are coverage-qualified
coordinate comparisons, not established mechanisms or independent evidence for
selection. Anchor eligibility alone leaves outside regions variably represented.

Full tables are in
`results/structural_comparisons/domain-anchored-coverage-20260927-v1`;
versioned completion evidence is in
`metadata/domain_anchored_coverage_completed_20260927.json`.

## Paired duplicate-to-reference contrasts

The [three-panel figure](figures/domain_anchored_contrasts_20260927.pdf) shows all
13 cases, with n50/c70 coverage qualification indicated by color. Axes have
separate scales. The intervals cover every mapping/mask/boundary alternative,
not statistical uncertainty. All 832 A-reference/B-reference pairs were joined
by their exact configuration **before** subtraction; an independent full CSV
reconstruction checked each contrast. All 234 case/screen/measurement rows and
all alternatives remain available.

For the three cases with both anchor and outside coverage passing n50/c70,
signed A-reference minus B-reference residual ranges (Å) are:

| Case | Domain core | Outside under domain transform | Outside fitted independently |
| --- | ---: | ---: | ---: |
| Heliocybe OG0000054 | +0.162 to +0.165 | −0.876 to −0.542 | −0.643 to −0.194 |
| Jaapia OG0000294 | +0.162 to +0.187 | −1.069 to −0.880 | −0.183 to +0.159 |
| Phycomyces OG0002812 | −0.130 to −0.129 | +0.577 to +0.603 | +0.109 to +0.122 |

The Jaapia contrast loses a stable direction when the outside coordinates are
fitted independently. This motivates testing whether relative arrangement
contributes to the apparent asymmetry; it does not identify a particular
interdomain rotation or establish its biological reality. The outside sets can
include flexible linkers and multiple domains, and prediction uncertainty remains
unquantified. Heliocybe and Phycomyces retain outside direction under either fit,
so the contrast cannot be attributed solely to fixing the domain transformation.

Run `scripts/summarize_domain_anchored_contrasts.py` to recreate tables and figures.
Outputs are in
`results/structural_comparisons/domain-anchored-contrast-summary-20260927-v1`;
completion and visual inspection evidence is in
`metadata/domain_anchored_contrast_summary_completed_20260927.json`.
