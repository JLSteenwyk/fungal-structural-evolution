# Combined domain coverage and geometry qualification

The queued join applies the same six original-interval coverage screens and
rotation-uniqueness criterion to both target and background domain cohorts.
It retains every pair under both masks, including unavailable inputs and all
coverage, RMSD and nonunique-rotation exclusions. Both alignment input orders
must pass. No structural metric is averaged or favorable order selected.

The inputs cover 70,395 target interval pairs and 66,929 background interval
pairs: 274,648 pair/mask rows in total. Target pairs include comparisons used
in event/reference analysis; these rows must still be projected to their
specific biological roles. They are not independent duplication events.
The join requires independently verified geometry and coverage receipts bound
to the same diagnostic source. All numeric-status and aligned-length bindings
must agree, and the entire geometry key set must be consumed exactly once.

The source plan is `metadata/domain_coverage_geometry_plan_20260927.json`;
output will be `results/structural_comparisons/domain-coverage-geometry-20260927-v1`.
The producer waits for the exact background geometry verifier. Independent
verification is queued separately and waits for successful producer termination.
It reconstructs every exported field using dataframe joins, compares complete
key sets, and checks all flags, exclusion strings and mask-specific/common
cohort counts. Production and full verification remain pending.

Each stage uses one CPU, 16 GiB RAM, no swap, and at most 2 GiB output; the
planning estimate is 1–15 minutes per stage after dependencies finish.
No GPU predictions or paid resources are used. Scripts parse successfully;
full-data verification is required before reporting the joined results.

Passing these criteria establishes numerical identifiability and the stated
coverage only. Domain boundary uncertainty, prediction uncertainty, matched
event eligibility, repeated control use, family and phylogenetic dependence
remain necessary for duplication-effect inference. The original failed strict
RMSD audits remain preserved and are not reclassified as passing.
