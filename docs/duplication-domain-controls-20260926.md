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

Production completed successfully. Its immutable receipt says
`pending_readback`; the separate full independent readback below now passed.
There are 542,964 pair/policy rows, with 37,755, 37,750, 37,595 and 37,590
candidate single-copy Domain matches under alignment E-value, alignment
bit-score, envelope E-value and envelope bit-score policies respectively.
These overlapping policy counts are not additive. Among the 227,089 models,
96,051 have no qualifying annotation and 1,944 have policy disagreement.
These producer counts now have a complete independent readback; all source
and output pins were checked again after readback completion.

## Independent complete output readback passed

`scripts/readback_duplication_domain_controls.py` now reconstructs every model
annotation export from separate registry segment and policy-membership queries,
then reconstructs the full primary/reference pair union, each architecture
category and every candidate single-copy domain pair. It does not import the
producer's classification or matching functions. The independent matcher uses
retained Domain occurrence counts and pairwise candidate matching, preserving
ineligible repeats when assessing copy uniqueness. All output fields, source
identities, complete model/pair/policy universes and aggregate counts must agree.

Known empty/single-copy/repeated-domain cases and 1,000 varied agreement cases
passed with `python scripts/check_duplication_domain_control_readback.py`.
The full run uses one CPU, 16 GiB RAM and no swap, GPUs or paid resources;
0.1–6 hours is an uncalibrated planning interval. Its plan and exact process
identity are in `metadata/duplication_domain_control_readback_plan_20260926.json`
and the corresponding launch record. Run with:

```bash
python scripts/readback_duplication_domain_controls.py \
  --plan metadata/duplication_domain_control_readback_plan_20260926.json
```

The output is
`results/structural_comparisons/duplication-domain-control-readback-20260926-v1.json`.
The readback passed all 227,089 models, 135,741 pairs, 542,964 policy rows
and all 150,690 domain-match policy records (overlapping across policies).
Its receipt is archived in
`metadata/duplication_domain_control_readback_completed_20260926.json`.
It checks interpretation of the same
frozen annotations, not an independent Pfam search or structural validation of
domain boundaries.

## Full domain-interval comparison workload

`scripts/prepare_duplication_domain_pairs.py` expanded every verified candidate
Domain match to both alignment and envelope boundaries. It retained all
301,380 policy/boundary links and deduplicated only identical model-coordinate
intervals and interval pairs. The resulting workload contains 114,590 intervals
from 50,204 models and 70,395 distinct interval pairs. Intervals contain
19,891,530 residues in total, with maximum length 1,066. Full and pLDDT>=70
comparisons in both input orders would require 281,580 dispositions.

The link table retains source protein pair, scope, policy, accession/version,
hit identifiers and boundary choice. Interval identities include model,
version, sequence hash and exact inclusive bounds. Identical alignment and
envelope bounds share computations while retaining both source links.

Plan: `metadata/duplication_domain_pair_plan_20260926.json`.
Launch: `metadata/duplication_domain_pair_launch_20260926.json`.
Completed producer record: `metadata/duplication_domain_pair_completed_20260926.json`.
Output: `results/structural_comparisons/duplication-domain-pair-inventory-20260926-v1/`.
The inventory used one CPU and 16 GiB RAM, no swap or GPUs. Its planning output
allowance was 1 GiB with uncalibrated 0.01–2 hour runtime. Reproduce with:

```bash
python scripts/prepare_duplication_domain_pairs.py \
  --plan metadata/duplication_domain_pair_plan_20260926.json
```

Interval identity/version/boundary distinction and out-of-range rejection were
checked before execution. All source and output hashes were checked after
completion. The producer status remains pending independent full link/interval
readback. No domain coordinates or alignments have been generated by this
inventory; confidence masks, PAE/orientation controls and biological inference
remain separate requirements.
