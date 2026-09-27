# Shared measurements across selected gene families

The full dependency inventory and independent verifier passed. Starting from
all 2,786,912 selected records, the inventory takes the union of participating
nodes across scenarios/policies, separately within each guide. It includes
19,250 target nodes per guide, 13,393 profile-guide backgrounds and 13,394
MAFFT-guide backgrounds. All 391,722 endpoint/entity rows were reconstructed
from native selected-node sources. This union precedes structural qualification.

For each guide, three graphs connect families sharing an exact gene ID, model
ID/version or sequence hash; a fourth combines those identity types. Both guides
contain 2,696 selected families. No gene identifier connects distinct families.
Model identities and sequence hashes each connect exactly one pair of families,
OG0000004 and OG0002812, leaving 2,695 combined components. The maximum component
contains two families; all other families are singletons under these identities.
The combined graph has four entity links across that family pair, representing
two models and their two sequence hashes, not four independent observations.

Shared models are AF-A0A2U1IYQ5-F1 and AF-A0A2U1J4X4-F1, version 6. Their endpoint
records involve genes F133377_PVZ97944.1, F133377_PWA00111.1,
F61424_PVU89359.1 and F61424_PVU93954.1 in taxa F133377 and F61424. Exact rows,
roles, guide assignments and sequence hashes are in
`metadata/selected_control_cross_family_endpoints_20260927.json`.
This is shared-measurement dependence, not proof that families should be merged
biologically or that an annotation is erroneous. Orthology assignments remain
unchanged. Downstream family-based uncertainty analysis must retain this link
rather than assume these two families have entirely independent measurements.

The producer uses union-find to form components. The independent checker
reconstructs every endpoint tuple from all selected-node sources, checks the
complete cross-family entity list, and rebuilds every component using SciPy
sparse graph traversal. Exact node/row counts, memberships and summaries agree.
Both services terminated successfully. Phylogenetic dependence between other
families/taxa is not tested or removed by this identity graph.

Plan: `metadata/selected_control_family_dependency_plan_20260927.json`.
Output: `results/orthology/selected-control-family-dependencies-20260927-v1`.
Proof: `metadata/selected_control_family_dependency_readback_20260927.json`.
Scripts: `audit_selected_control_family_dependencies.py` and
`readback_selected_control_family_dependencies.py`.
Each service used one CPU, 8 GiB RAM and no swap, with 0.1–2 hours planned and
2 GiB/1 GiB output allowances. No GPU predictions or paid resources.

All original gene-family labels and scenarios remain available. Identity-based
components are potential dependence blocks for sensitivity analysis, not a
validated complete covariance model or a claim of independent evolutionary
replicates. Structural qualification can remove records and reduce dependence;
this complete selected-node union preserves every possible shared-identity link
before that filtering.
