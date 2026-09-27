# Matched domain qualification grid

The full producer completed 3,777,048 configuration/boundary/mask/screen cells
for all 104,918 distinct matched protein-pair/policy configurations. Independent
full verification is running; the output remains pending acceptance until its
proof is complete. Source inventories contain 150,080 shared boundary-domain
matches, including repeated policies and pair configurations.

For each shared Pfam accession, require both the target and background domain
comparison to pass the same coverage and numerical-geometry screen. Each input
comparison already requires both alignment orders. Preserve alignment and
envelope boundaries separately, all six screens, and full, pLDDT70 and both-mask
cohorts. The both-mask cohort is the intersection of eligible Pfam accessions,
not merely configurations with some passing domain in each mask.

Each cell retains the original configuration status, shared-domain count,
eligible and ineligible accession lists, and a status for all/some/no shared
domains passing. Identical-model configurations remain explicitly deferred,
without invented zero scores. Configurations with no shared comparison remain
in every cell; this is not evidence that the biological domains are absent.
No favorable threshold, mask, domain or boundary is selected.

Plan: `metadata/selected_control_domain_qualification_plan_20260927.json`.
Output: `results/structural_comparisons/selected-control-domain-qualification-20260927-v1`.
Producer: `scripts/qualify_selected_control_domains.py`.
Verifier: `scripts/readback_selected_control_domain_qualification.py`.
The verifier independently intersects target/background pass sets and checks
all cells, exact accession lists, missingness categories and aggregate counts.
It waits for the exact producer identity and successful service termination.
Final proof will be `metadata/selected_control_domain_qualification_readback_20260927.json`.

Both stages use one CPU, 16 GiB RAM and no swap, with 0.1–2 hours budgeted per
stage. Output budgets are 4 GiB and 1 GiB. No GPU prediction or paid resources.
Source files and code are checksum-bound. These are repeated configuration
counts rather than independent evolutionary events; downstream projection must
retain all 2,786,912 selected records and account for shared controls and families.
Qualification is conditional on the observed models and chosen matching design.
It does not establish prediction accuracy, calibrated uncertainty or a causal
relationship between duplication and structural divergence.
