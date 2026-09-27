# Whole-protein common-residue comparisons

The full producer completed 563,808 fit records for 17,619 oriented protein-model
triads: two confidence masks, eight combinations of native alignment order and
two residue-mapping definitions. Each record compares duplicate A, duplicate B
and their reference using identical residue triples for all three distances.
Reference-intersection and cycle-consistent mappings remain separate.

Independent validation is running. It reconstructs every record from the
hashed coordinate inputs and residue maps, checks quaternion fits against the
producer's SVD fits, and checks confidence, sequence identity, geometry,
exclusions and all six coverage screens. Production completion alone does not
establish validation success.

The queued `scripts/summarize_verified_whole_protein_common_fits.py` waits for
the exact audit process and requires terminal success, matching source hashes
and all 563,808 checked records. It then groups the eight alignment orders for
each triad/mask/mapping definition. For each minimum length (30 or 50 residues)
and minimum original-protein coverage (50%, 70% or 90%), it counts triplets
passing all eight orders, some orders or no orders. It retains separate strata
for one, two and three distinct source models. A second complete CSV-based
aggregation must reproduce the pandas aggregation and the serialized summary.

The expected output is
`results/structural_comparisons/whole-protein-common-fit-summary-20260927-v1/order_robust_coverage.tsv`.
Completion requires
`metadata/whole_protein_common_fits_completed_20260927.json` and a successful
terminal summary service; neither is asserted by this launch-stage document.
The summary uses one CPU and at most 8 GiB RAM, with a planning allowance of
approximately 1–30 minutes after the audit. It does not start GPU inference.

These are model-triad counts, not independent duplication events. A model may
be reused across events or references, and shared-model triads require separate
interpretation. Coverage robustness does not establish biological asymmetry,
ancestral states or a significant duplication effect. Event/reference links,
phylogenetic dependence, prediction uncertainty and background comparisons
remain necessary for the evolutionary analysis.
