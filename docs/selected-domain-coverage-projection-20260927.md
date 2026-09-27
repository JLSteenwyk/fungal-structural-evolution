# Domain availability among selected matched records

The projection joins all 2,786,912 selected records to the fully verified
configuration qualification grid. It covers 100,328,832 record/boundary/mask/
screen combinations across 54 matching scenarios, four annotation policies,
two phylogenetic guides, two domain boundaries, three mask cohorts and six
coverage screens. All 15,552 summary strata must be present, including zeros.
Production and independent verification are running/queued, not yet complete.

Each stratum reports the metadata-selected denominator, usable selected records,
eligible domain occurrences, all five availability statuses, selected and usable
taxa/families/background nodes, and maximum reuse of a usable background node.
The denominator is selected controls; unmatched targets remain in the upstream
selection inventory and must also enter overall project coverage reporting.
Configuration counts are not substituted for selected-record counts.

The producer uses DuckDB SQL joins and grouped distinct counts. The independent
verifier reads the full selected records, target metadata and qualification
cells, builds compact NumPy arrays, and reconstructs all counts, distinct sets
and reuse maxima with pandas/NumPy without using the producer database or helpers.
It requires exact producer identity and successful termination, complete bound
source proofs, all input hashes, and an exact full-stratum match.

Plan: `metadata/selected_domain_coverage_projection_plan_20260927.json`.
Output: `results/structural_comparisons/selected-domain-coverage-projection-20260927-v1`.
Producer: `scripts/project_selected_domain_coverage.py`.
Verifier: `scripts/readback_selected_domain_coverage.py`.
Final proof: `metadata/selected_domain_coverage_projection_readback_20260927.json`.
The DuckDB working database is a derived intermediate, not a substitute for
checksum-bound exported TSV and source manifests.

Each service uses one CPU, 16 GiB RAM and no swap. The producer limits DuckDB
to 12 GiB memory and 20 GiB temporary disk; output allowance is 8 GiB, with
0.1–4 hours planned per stage after dependency completion. The verifier has a
1 GiB output allowance. No GPU predictions or paid infrastructure are enabled.

These are descriptive availability and dependence diagnostics. Guide, policy,
scenario, boundary and mask alternatives are dependent sensitivity analyses.
Usable counts do not establish independent replicates, calibrated prediction
uncertainty, selection neutrality or a causal effect of duplication. Structural
outcomes and biological effect models remain downstream.
