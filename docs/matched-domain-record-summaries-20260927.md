# Record-level structural contrasts and weighting sensitivity

The full summary stage and independent verification completed across 192 boundary/mask/cohort/threshold/
input-order combinations, each with all 432 guide/policy/matching-scenario
strata: 82,944 summary rows. No input order or sensitivity setting is
averaged with another. All 2,786,912 metadata-selected records remain represented
in denominators; unsupported strata have zero counts and blank means.

Within a record, first take an equal mean over its eligible shared Pfam domains.
This defines an explicit record-level estimand: proteins with more eligible
domains do not automatically contribute more weight to the record-weighted
summary. Exact per-configuration means and eligible-domain counts are stored
in 192 checksum-bound Parquet partitions, preserving the source domain links
for alternative future weighting. These are simple domain averages, not common
four-protein residue fits, inferred evolutionary rates or physical displacements.

Report nine descriptive measurements: target-minus-background RMSD and sequence
identity; each side's RMSD and identity; difference in minimum original-interval
coverage; log aligned-length ratio; and difference in joint high-confidence
residue fraction. For each measurement, report three separate views: equal
weight per matched record, equal weight per family after within-family record
averaging, and equal weight per taxon after within-taxon record averaging.
These views describe sampling sensitivity and are not phylogenetic corrections.
No p-values, confidence intervals or causal duplication effects are produced.

Every stratum's eligible-record, family, taxon and background counts must agree
exactly with the completed independently verified coverage projection. Full
independent verification is queued: it will rebuild every configuration mean
from all raw matched measurements and recompute each record/family/taxon summary
using pandas/NumPy without the SQL producer or helper functions. Empty means,
all counts, denominator bindings and the exact 192-setting grid are included.

Plan: `metadata/matched_domain_record_summary_plan_20260927.json`.
Output: `results/structural_comparisons/matched-domain-record-summaries-20260927-v1`.
Producer: `scripts/summarize_matched_domain_records.py`.
Verifier: `scripts/readback_matched_domain_record_summaries.py`.
Final proof: `metadata/matched_domain_record_summary_readback_20260927.json`.
Production and verification passed in full: 99,830,400 configuration values
were independently reconstructed across all 192 settings, and all 82,944
record/family/taxon summaries matched. Every configuration-partition hash was
rechecked after both services terminated successfully. Completion is archived in
`metadata/matched_domain_record_summary_completed_20260927.json`.

Each stage reserves one CPU, 16 GiB RAM, no swap, and 1–8 hours after dependencies.
The producer limits DuckDB to 12 GiB memory and 20 GiB temporary disk and budgets
16 GiB output; the verifier budgets 1 GiB output. Every setting is processed
separately to bound aggregation memory. No predictions or paid resources.

After verification, inspect sequence/coverage imbalance and weighting/order
sensitivity before fitting family- and phylogeny-aware models. Shared controls,
shared ancestry, prediction uncertainty, conditional cohort selection and the
fact that separate alignments can map different residues remain unresolved by
these descriptive means.


For a traceable illustrative cell (profile guide, alignment-evalue policy,
alignment boundaries, pLDDT70 measurements restricted to the both-mask cohort,
n50/c70, input orders 0/0), S45 has 6,301 matched records. Its raw target-minus-
background RMSD mean is 0.040944 Å with record weighting, 0.018241 Å with equal
family weighting and 0.092330 Å with equal taxon weighting. Its mean sequence
identity difference is -0.005873 and original-coverage difference -0.008670.
The corresponding focal-only S46 cell has 495 records and RMSD means 0.063658,
0.033379 and 0.042936 Å; identity difference is -0.011501 and coverage difference
-0.007302. These rounded examples are descriptive and unadjusted, not evidence
of statistical significance or a causal duplication effect. The remaining
orders and complete grid are retained in the source and sensitivity report.
