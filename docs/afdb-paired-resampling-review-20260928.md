# Expanded AlphaFold paired resampling

The 125-marker paired AA/3Di resampling and its native-output audit completed
successfully. The controller is inactive with exit status zero. Its receipt
binds the audited outputs to the versioned production plan; both output-table
hashes were checked again on September 28.

- 375 marker/block batches: 200 draws each at block lengths 1, 10 and 30.
- 75,000 attempted paired draws; 74,977 estimable and 23 unestimable.
- 149,954 native fitted trees validated by the producer's audit.
- 371,358 conditional interval rows, covering 61,893 branches per fit type
  and block setting.
- Every batch retains at least 197 of its 200 draws. The 23 unestimable draws
  occur in 13 batches; all rows meet the predeclared 90% estimability threshold.

The completed audit checks deterministic resampling columns and FASTAs,
source/model/topology bindings, native tree branch lengths and all recorded
artifact hashes. A separate September 28 review additionally checks table
uniqueness, row counts, quantile ordering, finite values and paired covariance
consistency, and inventories warnings directly from receipt-bound native
reports and logs. That warning review completed successfully, with all 149,954 fits inventoried.

Warnings occur in 148,752 of the 149,954 native fits. The complete census found:

| Warning | Fits containing warning |
| --- | ---: |
| Near-zero internal branches | 148,324 |
| Rare alignment states | 47,505 |
| Sequences with over 50% gaps/ambiguity | 90,732 |
| Very long branches (>9.8 substitutions/site) | 912 |

Classes overlap. All warning lines are classified; classification does not
resolve their scientific implications. All 23 unestimable draws contain an
entirely missing taxon after resampling. No failed draw was silently replaced.

Across block settings, 8,930–9,099 of 61,893 structural-alphabet branches and
5,023–5,050 amino-acid branches have at least half their resampled estimates
at or below 1e-5 substitutions per site. These counts include terminal and
internal branches and are separate from IQ-TREE's internal-branch warning
threshold. They describe limited branch information; they do not establish
absence of structural evolution or define accelerated branches.

Percentile intervals remain conditional on the fixed sequence topology and
specified models. They do not incorporate all topology, model, prediction or
spatial-dependence uncertainty. Paired sampling covariance concerns estimation
error, not evolutionary sequence–structure coupling. No ratios, significance
claims or acceleration rankings follow from this checkpoint alone.

Evidence and reproduction:

- `metadata/recovered_afdb_paired_resampling_plan_20260926.json`
- `metadata/recovered_afdb_paired_resampling_completed_20260928.json`
- `metadata/afdb_paired_warning_review_launch_20260928.json`
- `metadata/afdb_paired_warning_review_completed_20260928.json`
- `metadata/afdb_paired_resampling_warning_types_20260928.json`
- `results/phylogeny/paired-resampling-audit-afdb-recovered-20260926-v1`
- `scripts/review_afdb_paired_resampling_outputs.py`

The separate review writes to
`results/phylogeny/paired-resampling-review-afdb-recovered-20260928-v2`.
An incomplete first review was intentionally terminated to correct its
warning-classification pattern; its outputs remain preserved and are not
used as completed evidence. The underlying resampling was not restarted.
