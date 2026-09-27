# Whole-protein common-residue comparisons

The full producer completed 563,808 fit records for 17,619 oriented protein-model
triads: two confidence masks, eight combinations of native alignment order and
two residue-mapping definitions. Each record compares duplicate A, duplicate B
and their reference using identical residue triples for all three distances.
Reference-intersection and cycle-consistent mappings remain separate.

Independent validation passed all 563,808 records reconstructed from hashed
coordinate inputs and residue maps. Quaternion and SVD RMSDs/contrasts agreed
within 2.42 × 10⁻¹³ Å. All confidence, sequence identity, geometry, exclusion
and six coverage-screen values passed the full readback.

`scripts/summarize_verified_whole_protein_common_fits.py` completed after the
exact audit process terminated successfully. Both full-table aggregations and
serialized readback agreed across 70,476 triad/mask/mapping groups and 72
coverage-summary rows. Results retain separate one-, two- and three-model
strata, all eight alignment orders, and six length/coverage thresholds.

`scripts/link_whole_protein_common_fit_eligibility.py` then linked all 36,944
event/reference records without selecting among tied references. The table
retains all triplets, including shared-model and failing cases. The following
counts require three distinct models and a passing reference in **all 32**
mask/order/mapping configurations:

| Minimum residues / original-protein coverage | MAFFT-guide events | Profile-guide events |
| --- | ---: | ---: |
| 30 / 50% | 5,877 | 5,878 |
| 30 / 70% | 3,538 | 3,538 |
| 30 / 90% | 948 | 948 |
| 50 / 50% | 5,851 | 5,852 |
| 50 / 70% | 3,534 | 3,534 |
| 50 / 90% | 948 | 948 |

Denominators are 17,448 MAFFT-guide events and 17,461 profile-guide events in
this modeled-reference inventory. Counts are per guide and must not be added
across guides. At least one reference passing does not imply agreement across
all tied references. Matching aggregate counts do not establish identical
event membership across guides.

Outputs reside in `results/structural_comparisons/` under
`whole-protein-common-fit-summary-20260927-v1` and
`whole-protein-common-fit-event-links-20260927-v1`. The latter includes full
triad eligibility, every event/reference link and the guide/threshold summary.
Hashes and validation evidence are recorded in
`metadata/whole_protein_common_fits_completed_20260927.json` and
`metadata/whole_protein_common_fit_event_links_completed_20260927.json`.
The analyses used one CPU each and did not start GPU inference.

These are model-triad counts, not independent duplication events. A model may
be reused across events or references, and shared-model triads require separate
interpretation. Coverage robustness does not establish biological asymmetry,
ancestral states or a significant duplication effect. Event/reference links,
phylogenetic dependence, prediction uncertainty and background comparisons
remain necessary for the evolutionary analysis.
