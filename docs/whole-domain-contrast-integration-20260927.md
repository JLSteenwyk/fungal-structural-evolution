# Whole-protein and domain contrast integration

All 77,760 domain-comparison rows were joined to exact family/gene-pair/screen
whole-protein comparisons. The integration retains six length/coverage screens
and three descriptive margins (0, 0.01 and 0.1 Å). Whole-protein directions were
recomputed at each domain margin; the earlier numerical-zero classification
was not substituted for these margins.

For each guide, reference-gene sets and all A/B/reference model accession/version
combinations must match exactly before the comparison receives an agreement or
opposition label. Whole-protein comparisons must have every reference eligible;
domain comparisons retain their original requirement for all references,
annotation policies, boundaries and fit settings. Both guides must support the
same direction within each scale. Missingness and reference/model-set mismatches
remain explicit.

At the 50-residue/70%-coverage screen and the **descriptive 0.1 Å margin**:

| Relationship | Event/domain combinations |
| --- | ---: |
| Same structural direction at both scales | 320 |
| Opposite structural directions between scales | 6 |
| Domain direction incomplete or unresolved | 302 |
| Whole-protein direction incomplete or unresolved | 3,516 |
| Reference or model sets differ | 176 |

These 4,320 combinations are the domain-comparison inventory for this screen,
not independent events. The whole-protein universe additionally contains
100,014 pair/screen rows with no domain-comparison row across the six screens;
these are exported separately and must not be interpreted as domain absence.
All margins and screens remain in the output; the table above is one stated
view of the full result.

Opposing directions identify comparisons for closer inspection. They do not
by themselves establish altered domain orientation: the residue sets, fitting
objectives and coverage denominators differ between whole-protein and domain
comparisons. Domain composition, flexible linkers, confidence, annotation and
prediction-source effects also require review. Whole-protein coverage is measured
against original protein lengths; domain coverage is measured against domain
lengths. The analysis does not subtract RMSDs across these scales or claim a
statistical duplication effect.

Run `scripts/integrate_whole_and_domain_contrasts.py` using one CPU. It checks
source hashes, independently recomputes all whole-protein direction labels and
reads back every serialized output cell. The domain source has a full prior
independent readback. Outputs and receipt are in
`results/structural_comparisons/whole-domain-contrast-integration-20260927-v1`;
versioned provenance is in
`metadata/whole_domain_contrast_integration_completed_20260927.json`.
