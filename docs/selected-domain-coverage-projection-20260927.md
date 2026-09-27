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

Current plan: `metadata/selected_domain_coverage_projection_plan_20260927_v2.json`.
Current output: `results/structural_comparisons/selected-domain-coverage-projection-20260927-v2`.
Current producer: `scripts/project_selected_domain_coverage_v2.py`.
Verifier: `scripts/readback_selected_domain_coverage.py`.
Final proof: `metadata/selected_domain_coverage_projection_readback_20260927_v2.json`.
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


## Full-scope resource recovery

The first producer terminated with exit code 1 after its full-grid grouped
join reached DuckDB's explicit 20 GB temporary-disk allowance. Its checker
also terminated with exit code 1 at the unsuccessful-producer gate. No
completed receipt was produced; original scripts, plan, output database and
failed services are retained. The failure and terminal states are archived in
`metadata/selected_domain_projection_failure_recovery_20260927.json`.

The v2 producer uses identical aggregate expressions but processes one
boundary/mask/screen setting at a time. All 36 settings, 100,328,832 selected
record cells and 15,552 summary strata remain in scope. Resource limits are
unchanged. The unchanged independent verifier is queued with the new plan,
exact process identity and a new output path. Fourteen settings completed at
the live checkpoint, with cgroup memory near 0.8 GiB. This is evidence of
progress, not completion or a full-memory benchmark.
