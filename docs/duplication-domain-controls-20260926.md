# Domain controls for duplication and reference comparisons

The full primary/reference union comprises 227,089 models and 135,741 distinct
model pairs. The completed inventory joins all models by exact sequence,
model identifier/version, length and path to the independently audited
whole-proteome structural domain registry. It preserves all four annotation
competition policies and all retained Pfam hit occurrences and types.

`model_annotations.jsonl` retains the source interval identities, boundaries,
HMM coverage and conservative/candidate flags for each model and policy.
`pair_architecture_controls.tsv` compares ordered accession/version/type
signatures and distinguishes:

- Neither model annotated.
- One model unannotated.
- Same ordered annotations, including copy multiplicity.
- Same annotation content and multiplicity but different order.
- Different annotation content or multiplicity.

Rows also retain each model's conservative-architecture status, Domain hit
counts and candidate single-copy domain-match counts. Unannotated does not
mean biologically absent. Order/content differences are annotation contrasts,
not inferred gain, loss, fusion or rearrangement events. Same-order proteins
can still differ in domain orientation or unannotated regions.

`single_copy_domain_pairs.tsv` identifies exact Pfam-accession/version matches
where both proteins have exactly one retained occurrence of that Domain and
both occurrences have the registry's conservative candidate flag. Counting
all retained occurrences prevents a partial/ineligible second copy from making
a repeated domain appear unique. Pfam Family and other non-Domain types remain
in the architecture signatures but are not promoted to Domain intervals.
Matches retain original source-protein alignment boundaries. They still need
coordinate extraction, residue confidence and orientation/PAE assessment,
structural alignment and biological interpretation. Policy rows overlap and
must not be summed as independent comparisons.

Script: `scripts/inventory_duplication_domain_controls.py`.
Plan: `metadata/duplication_domain_control_plan_20260926.json`.
Completed run record: `metadata/duplication_domain_control_completed_20260926.json`.
The job finished before live process-identity capture; the record includes
the journal PID, successful terminal service state and checked receipt hashes.
Output: `results/structural_comparisons/duplication-domain-controls-20260926-v1/`.

Resources are one CPU, 16 GiB RAM, no swap and an estimated 2 GiB output. The
0.1–6 hour planning interval is uncalibrated. No GPU or paid services are used.
Reproduce with:

```bash
python scripts/inventory_duplication_domain_controls.py \
  --plan metadata/duplication_domain_control_plan_20260926.json
```

`python scripts/check_duplication_domain_controls.py` passed all five signature
categories, preservation of repeated/versioned annotations, exclusion of
non-Domain and ineligible hits, and prevention of false single-copy matching
when an ineligible additional domain copy is retained.

Production completed successfully. The receipt explicitly says
`pending_readback`; independent full output verification remains required.
There are 542,964 pair/policy rows, with 37,755, 37,750, 37,595 and 37,590
candidate single-copy Domain matches under alignment E-value, alignment
bit-score, envelope E-value and envelope bit-score policies respectively.
These overlapping policy counts are not additive. Among the 227,089 models,
96,051 have no qualifying annotation and 1,944 have policy disagreement.
These are producer counts pending independent readback; all source pins and
output hashes were checked after successful completion.
