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
completion. The immutable producer status says pending readback; the separate full
link/interval readback below has now passed. No domain coordinates or alignments have been generated by this
inventory; confidence masks, PAE/orientation controls and biological inference
remain separate requirements.

## Complete interval-workload readback passed

`scripts/readback_duplication_domain_pairs.py` independently reconstructs all
301,380 source-policy/boundary links, checks every original source-match field,
and verifies both endpoints against the corresponding alignment/envelope
bounds in the audited annotations. It also checks all 114,590 interval
provenance records against source model descriptors and verifies the exact
70,395-pair deduplicated set, interval-use set, identities and summary counts.
No producer interval/pair construction function is imported.

The run uses one CPU, 16 GiB RAM, no swap, GPU or paid resources, with an
uncalibrated 0.01–2 hour planning interval. The plan is
`metadata/duplication_domain_pair_readback_plan_20260926.json`. The process
completed before live identity capture; the completed archive records the
journal PID and receipt proof. Reproduce with:

```bash
python scripts/readback_duplication_domain_pairs.py \
  --plan metadata/duplication_domain_pair_readback_plan_20260926.json
```

Output is
`results/structural_comparisons/duplication-domain-pair-readback-20260926-v1.json`.
The full readback passed all links, intervals, pairs and aggregate totals.
Its receipt and post-completion pin verification are archived in
`metadata/duplication_domain_pair_readback_completed_20260926.json`. This verifies
interval bookkeeping against the same annotations, not physical domain
boundaries or biological evolutionary events.

## Full domain coordinate preparation queued

`scripts/materialize_duplication_domain_inputs.py` waits for both exact
full-model coordinate-readback processes and requires their successful bound
receipts. Every source shard must match its hash and independent proof. It
then extracts all 114,590 verified intervals from 50,204 models, producing
229,180 full/pLDDT>=70 dispositions. Source path, hash, sequence hash and model
length are checked for every interval. Short masks and source rejections remain
explicit. All expected intervals must appear exactly once per mask.

The renderer in `scripts/duplication_domain_inputs.py` preserves original
full-protein residue numbers through interval slicing and masking. It uses
the same checked C-alpha serializer as the whole-protein pipeline, with XYZ
rounded to 0.001 Å and confidence to 0.01. Independent Bio.PDB parsing fixtures
checked original numbering, coordinates and confidence for full and sparse
masks; empty masks and invalid bounds were checked separately. A completed
handoff fixture covered both source collections, rejected models, written
hashes and altered-proof rejection. Reproduce these checks with:

```bash
python scripts/check_duplication_domain_inputs.py
```

Plan: `metadata/duplication_domain_input_plan_20260926.json`.
Exact launch identity: `metadata/duplication_domain_input_launch_20260926.json`.
Output: `results/structural_comparisons/duplication-domain-inputs-20260926-v1/`.
Run with `python scripts/materialize_duplication_domain_inputs.py --plan
metadata/duplication_domain_input_plan_20260926.json`.

The queued job uses one CPU, 8 GiB RAM, no swap, an 8 GiB output planning
allowance and 100 GiB free-disk reserve. Its uncalibrated planning interval is
0.1–8 hours after dependencies finish. No GPU or paid resources are used.
Production extraction is waiting, not completed. Independent full domain
serialization readback and domain alignments remain to be prepared. No PAE,
biological boundary or evolutionary-event validation is implied.

## Independent domain-coordinate readback queued

`scripts/readback_duplication_domain_inputs.py` now waits for the exact domain
extraction process and requires its successful plan-bound receipt. It checks
the complete 229,180 interval/mask grid, source model/shard identities, upstream
coordinate receipt/proof bindings and all ready PDB file hashes. It independently
selects residues from the audited full-model arrays and checks every PDB atom's
original residue number, amino-acid identity, XYZ, confidence, occupancy and
serialization fields. It does not import the extraction renderer. Coordinate
and confidence tolerances are fixed at their printed rounding limits:
0.000501 Å and 0.005001 confidence units. Short masks and rejected sources are
checked explicitly. Raw-CIF reconstruction is provided by the earlier source
coordinate audits, not repeated here.

The completed two-source extraction/readback fixture passed four dispositions
and ten PDB C-alpha atoms. A changed coordinate was rejected even after its
recorded PDB hash was updated. Reproduce with:

```bash
python scripts/check_duplication_domain_input_readback.py
```

Plan and exact process identity:
`metadata/duplication_domain_input_readback_plan_20260926.json` and the
corresponding launch record. Run with
`python scripts/readback_duplication_domain_inputs.py --plan
metadata/duplication_domain_input_readback_plan_20260926.json`.
Output will be
`results/structural_comparisons/duplication-domain-input-readback-20260926-v1.json`.

The queued readback uses one CPU, 8 GiB RAM and no swap, GPU or paid resources.
Its output is a small proof receipt; uncalibrated runtime planning is 0.1–8
hours after extraction. Production readback is waiting, not complete. Domain
alignment orchestration and biological analyses remain pending.

## Complete domain alignment workload queued

`scripts/run_duplication_domain_alignments.py` now waits for the exact full
domain-coordinate readback process. The handoff requires its passed receipt,
the matching materialization receipt/plan, the exact interval/mask grid and
all interval identities. The interval-pair table and canonical pair hashes
are checked before execution. Every checkpoint retains both interval IDs,
source model/version, bounds, mask, input file hashes, command, native output
and timing. The complete input manifest, producer receipt and passed readback
receipt jointly bind every checkpoint.

The scope is all 70,395 interval pairs × two orders × two masks, or 281,580
explicit dispositions. Alignment/envelope boundaries and all policy links
remain in the frozen inventory for later sensitivity analysis. Short masks,
rejected sources, native errors, parse errors and timeouts remain explicit;
no automatic substitution or retry occurs.

Resources: four CPU workers, 16 GiB RAM, no swap, 32 GiB output planning
allowance, 100 GiB free-disk reserve and 600-second per-call timeout. At most
64 futures are pending. The 8–480 hour planning range is uncalibrated and not
a guaranteed bound. No GPUs or paid resources are used. The native fixture
passed both input orders, masked exclusions, checkpoint reuse and hash checks,
plus rejection of a mismatched independent-input-readback receipt:

```bash
python scripts/check_duplication_domain_alignment_handoff.py
```

Plan: `metadata/duplication_domain_alignment_plan_20260926.json`.
Exact launch: `metadata/duplication_domain_alignment_launch_20260926.json`.
Output: `results/structural_comparisons/duplication-domain-alignments-20260926-v1/`.
Run with `python scripts/run_duplication_domain_alignments.py --plan
metadata/duplication_domain_alignment_plan_20260926.json`.

Production alignments are waiting for verified input preparation. The future
completion status explicitly retains `pending_readback`; independent alignment
numeric readback, boundary/annotation-policy sensitivity and biological tests
remain separate requirements.

## Full domain-alignment numeric readback queued

`scripts/readback_duplication_domain_alignments.py` waits for the exact domain
alignment process and its matching completed receipt. It reconstructs the
complete 281,580 interval-pair/order/mask grid from the frozen pair table,
checks interval/model/bound identities, checkpoint paths/hashes, native commands
and the complete audited-input bundle. No missing or repeated disposition is
allowed. Failed/excluded outcomes are checked for internal consistency rather
than rerun or independently adjudicated.

Every successful alignment is independently parsed and mapped to hashed PDB
residues; sequence identity, coverage and least-squares RMSD are reconstructed
with the numeric helper used by the whole-protein readbacks. TM-scores are
checked against native text only, not independently reoptimized. Original
protein residue numbering remains explicit. Confidence summaries use rounded
PDB values; no PAE or biological interpretation is supplied by the audit.

The native full-handoff fixture passed both orders and excluded masks; dropping
one excluded result was rejected despite an updated checkpoint-manifest hash:

```bash
python scripts/check_duplication_domain_alignment_readback.py
```

Plan and exact launch identity are in
`metadata/duplication_domain_alignment_readback_plan_20260926.json` and the
corresponding launch record. Output will be
`results/structural_comparisons/duplication-domain-alignment-readback-20260926-v1/`.
Run with `python scripts/readback_duplication_domain_alignments.py --plan
metadata/duplication_domain_alignment_readback_plan_20260926.json`.

Resources are one CPU, 16 GiB RAM, no swap, single-thread BLAS, a 256-input PDB
cache and 1 GiB estimated output. Uncalibrated runtime planning is 0.5–24 hours
after the producer. No GPUs or paid resources are used. The readback is queued,
not completed; policy/boundary sensitivity and biological tests remain pending.


## Complete duplicate/reference architecture controls

The complete provisional-reference ledger now has triad-level architecture
controls: 36,944 event/reference combinations across both guides, retaining
all nearest-reference ties. Each has four annotation-policy alternatives,
for 147,776 rows. Pairwise controls are joined for duplicate A/B, A/reference
and B/reference, preserving pair keys, model versions, ordered-annotation
classes and conservative flags. Original source ledgers and successful
independent readbacks are hash-bound before the join.

A triad's summary class follows this explicit priority: any identical-model
pair; otherwise any incomplete annotation; otherwise any content/order
difference; otherwise identical ordered annotations with at least one
nonconservative pair; otherwise conservative identical ordered annotations.
Every pair's underlying fields remain available even when a higher-priority
class determines the summary. Identity is not classified as missing annotation;
its conservative flag is blank because no distinct-pair control is consulted.

| Class, alignment/e-value policy | Profile | MAFFT |
| --- | ---: | ---: |
| Conservative matching ordered annotations | 5,169 | 5,164 |
| Matching ordered annotations, not conservative | 1,848 | 1,847 |
| Annotation content/order difference | 2,752 | 2,746 |
| Incomplete annotation | 6,433 | 6,436 |
| Identical model in triad | 2,276 | 2,273 |
| Total event/reference combinations | 18,478 | 18,466 |

Across all four policies, conservative matching counts range from 5,156–5,169
for profile and 5,151–5,164 for MAFFT. Of the 36,944 combinations, 332 change
summary class between policies and 36,612 retain the same class. These counts
include dependent guides and reference ties: they are not independent events
or a sample size for inference. Conservative matching Pfam annotations do not
establish domain homology, stable domain orientation, structural similarity,
reference orthology or an evolutionary gain/loss. Missing annotation is not
biological absence. These controls will stratify subsequent structural
comparisons; they are not tests of duplication-associated structural change.

Output: `results/structural_comparisons/duplication-triad-architecture-20260926-v1/`.
Independent dataframe joins checked the entire triad/policy grid, source
model and pair identities, all pair fields, identity dispositions, categories
and totals. Both receipts are archived in
`metadata/duplication_triad_architecture_completed_20260926.json`.
This table transformation uses one CPU, no new annotation search or GPU.
Reproduce in fresh output paths:

```bash
python scripts/prepare_duplication_triad_architecture.py \
  --controls results/structural_comparisons/duplication-domain-controls-20260926-v1 \
  --control-readback results/structural_comparisons/duplication-domain-control-readback-20260926-v1.json \
  --references results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1 \
  --reference-readback results/structural_comparisons/duplication-reference-comparison-readback-20260926-v1.json \
  --queue results/structural_comparisons/duplication-model-pair-queue-20260926-v1 \
  --output <fresh-output-directory>
python scripts/readback_duplication_triad_architecture.py \
  --source <fresh-output-directory> --output <fresh-readback.json>
```


## Full expanded primary duplication annotations (September 30)

The expanded September 28 queue contains 276,682 model/version identities and
134,812 distinct model pairs. All models, including endpoints of identical-model
events, were rejoined to the complete expanded domain registry under all four
original policies. The independent checker separately joins native segment and
policy-membership records, then reconstructs every pair category and candidate
single-copy domain match. All **539,248 pair/policy records** passed. Exact
source checks and completed producer/checker process journals are bound in
`metadata/expanded_duplication_domain_controls_completed_20260930.json`.

There are 126,713 models without retained Pfam hits and 2,353 with policy
disagreement. Under the alignment/e-value policy, 53,390 pairs have the same
ordered annotations, 13,411 have different annotation content, 20 have the same
content in a different order, 11,845 have one unannotated endpoint and 56,146
have neither endpoint annotated. Matching annotation order alone is not the
conservative-architecture criterion used for background eligibility.
Single-copy candidate matches number 32,347/32,329/32,202/32,184 under alignment
E-value/alignment bit score/envelope E-value/envelope bit score, respectively;
these overlapping policies are not independent observations.

Versioned `inventory_duplication_domain_controls_v2.py` and
`readback_duplication_domain_controls_v2.py` retain the original ordered,
repeat-aware annotation logic while reading the complete audited primary queue.
They require the exact producer/readback binding for all queue files. The older
producer's reference-specific additional-model partition is not reused against
a different base queue. This primary inventory does not refresh sibling-reference
comparisons or ancestral domain events; those still require their own full
expanded source ledgers. The expanded background/control workflow is documented
in [the refresh workflow](terminal-sister-backgrounds-20260927.md#expanded-catalog-refresh-september-30).
