# Structural coverage of terminal duplication candidates

`scripts/summarize_duplication_sampling_coverage.py` independently reconstructs
all reported Terminal events with exactly one gene on each side from both
complete native duplication tables. It joins both genes to the same frozen
AlphaFold structural bridge used by the structural-comparison queues and
compares event identity, singleton status and model-presence counts against
the original coverage exports. This includes events with no models or only one
model, which do not enter the current structural pair queue.

The output contains an event-disposition table and summaries by taxon and
family. Every sampled taxon is retained, including zero qualifying events.
Zero denominators produce blank fractions, not zero structural coverage.
Counts distinguish neither/one/both genes modeled and identical-model cases.
Species names, study roles and lineage strings come from the pinned sampling
manifest. Profiles and MAFFT are kept separate, not treated as independent
replicates. Reported terminal singleton-side events are a restricted candidate
class, not the full biological duplication history.

The immediate purpose is to quantify ascertainment of the modeled candidate
set before downstream divergence/asymmetry interpretation. Missing models are
not biological absences. Counts do not establish missing-at-random sampling,
correct ascertainment bias or account for phylogenetic dependence. The frozen
bridge excludes later retrievals and ESMFold structures; the table does not
claim current total structure availability. It also does not independently
validate every missing-model candidate's gene-tree node.

Plan: `metadata/duplication_sampling_coverage_plan_20260926.json`.
Launch: `metadata/duplication_sampling_coverage_launch_20260926.json`.
Output: `results/orthology/duplication-sampling-coverage-20260926-v1/`.
Run with:

```bash
python scripts/summarize_duplication_sampling_coverage.py \
  --plan metadata/duplication_sampling_coverage_plan_20260926.json
```

The job is running. Resources are one CPU, 8 GiB RAM, no swap and an estimated
1 GiB output. The 0.1–4 hour planning range is uncalibrated. No GPU or paid
resources are used. Completed results and independent aggregate readback remain
pending.
