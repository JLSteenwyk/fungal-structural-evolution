# Duplication-event structure coverage

A full inventory was launched September 26 to connect reported duplication
events to the existing frozen AlphaFold structure bridge. Both audited guides
are included: 1,204,638 profile events and 1,204,919 MAFFT events. All events,
including nested events, multifurcations with an empty second side and events
with no structures, receive an output row. Counts of genes across events
repeat descendants and are not independent observations.

The bridge contains 1,319,513 protein links. This inventory describes that
frozen source; it does not include subsequently retrieved models or the
separate ESMFold marker collection. Model availability has not yet been
qualified by confidence, domain coverage or alternative conformations.

The output records gene counts and model-covered gene counts for each event
side, while retaining family, species-tree node, gene-tree node, native
support and event type. A separate candidate table contains terminal events
with exactly one protein on each side and models for both. It retains exact
gene/model IDs and flags identical sequences and identical model IDs. No
native-support threshold is silently imposed; support is retained for later
review. Larger and incompletely covered events remain in the full inventory.

These candidates are for subsequent tree-membership validation, confidence
and alignment assessment, direct structural comparison and matched controls.
They do not by themselves establish a duplication-associated structural
effect. Asymmetry needs an appropriate ancestral/outgroup reference and
uncertainty treatment. Gene-tree membership and reconciliation uncertainty
remain requirements for event interpretation.

Implementation: `scripts/inventory_duplication_structure_coverage.py`. It
requires the successful full event-identity audit and independently validated
structure bridge, pins the exact source files, checks all hashes before/after
execution, and verifies event and nested gene-assignment totals against the
prior audit. Test with `python scripts/check_duplication_structure_coverage.py`: 
200 independent coverage fixtures passed and three malformed events were
rejected, including repeated genes. Tests cover missing models, complex and
empty-second-side events as well as simple terminal candidates.

Plan: `metadata/duplication_structure_coverage_plan_20260926.json`. Launch
identity: `metadata/duplication_structure_coverage_launch_20260926.json`.
Service: `fungal-duplication-structure-coverage-20260926.service`.
Output: `results/orthology/duplication-structure-coverage-20260926-v1/`.
Resources: one CPU, 8 GiB RAM, no swap, 2 GiB output allowance; an uncalibrated
0.1–4-hour planning range on existing local resources. No GPU work is involved.
The full inventory is running; successful production output and independent
readback remain pending.


### September 26: coverage production completed

The service exited successfully and produced the full inventory for both
guides. Source-plan binding and all artifact hashes were rechecked; receipt
archived in `metadata/duplication_structure_coverage_completed_20260926.json`.
The profile inventory contains 109,245 simple terminal two-model candidates
from 20,492 families; MAFFT contains 109,228 from 20,474 families. Each covers
153 taxa and flags 6,354 candidates with identical sequences/models. These
are guide-specific candidate counts, not independent or validated duplication
effects. Independent source-to-output readback, candidate overlap and
gene-tree membership validation remain pending.


## Independent candidate reconstruction and reported-node review

A full review of all 109,245 profile and 109,228 MAFFT candidates is now
running. `scripts/validate_duplication_structure_candidates.py` does not
import the coverage producer. It independently enumerates eligible terminal
singleton-side events from both complete native tables and the frozen model
database, checking every candidate field and rejecting omitted or extra
eligible candidates. This validates the candidate export, not every row of
the separate full event-coverage inventory.

For each candidate-bearing family, the review parses the previously audited
resolved gene tree and looks up the exact reported gene-node name. It records
whether that node has exactly the two expected protein descendants, and
separately whether they are its two direct tip children. Missing nodes,
missing families and descendant mismatches remain explicit output statuses;
a completed review is not automatically a passing result. Guide overlap uses
canonical pairs of full protein labels, without assuming that guide-specific
family names identify the same membership.

This establishes correspondence between reported candidate events and the
saved resolved trees. It does not independently infer duplication, validate
rooting/support, date events or establish structural divergence. Identical
sequence/model flags remain in output. Larger events remain in the full
inventory for later work.

Six tree-node fixtures passed, including swapped expected pair order, extra
descendants, absent nodes, unary branches, duplicate labels and missing node
names. Reproduce with `python scripts/check_duplication_candidate_tree_nodes.py`.
Run the full review with `python scripts/validate_duplication_structure_candidates.py
--plan metadata/duplication_candidate_tree_review_plan_20260926.json`.
The plan pins source event tables, model database, coverage receipt/export,
resolved trees and their successful full-membership readbacks. All pins are
checked before and after execution.

Service: `fungal-duplication-candidate-tree-review-20260926.service`.
Launch identity: `metadata/duplication_candidate_tree_review_launch_20260926.json`.
Output: `results/orthology/duplication-candidate-tree-review-20260926-v1/`.
Resources: one CPU, 8 GiB RAM, no swap and 1 GiB output allowance. The
0.1–6-hour planning interval is uncalibrated. Results remain pending.


## Completed candidate reconstruction and tree correspondence

The full candidate review completed successfully. Independent reconstruction
from native event tables and frozen models matched all 109,245 profile and
109,228 MAFFT candidate rows. Every reported node had exactly the two
expected proteins as its two direct tip children: no missing nodes, missing
family trees or descendant mismatches were found in this candidate set.

The guides share 108,918 exact protein pairs; 327 pairs are profile-only and
310 MAFFT-only. All shared pairs pass the node check in both guides. This
comparison uses protein identities, not family-name matching. Receipt,
source-plan binding and all output hashes were checked after service exit
zero and archived in
`metadata/duplication_candidate_tree_review_completed_20260926.json`.

These checks validate export eligibility and correspondence to saved trees,
not independent evidence for duplication biology, calibrated node support,
rooting, timing, structure quality or a structural effect. The complete
event-coverage inventory still needs separate full-row readback; only its
candidate export was independently reconstructed here.


## Direct structural-comparison queue

The complete reviewed candidate set has been joined to frozen model/version
records. The queue retains all 218,473 guide-specific event links; 12,708
rows reference the same model on both sides and remain explicitly flagged.
The remaining 205,765 event links map to 103,200 unique distinct model pairs,
requiring 206,400 alignments to evaluate both input orders. Shared model
comparisons can be computed once while retaining every gene, taxon, family
and guide association. This computational reuse does not make event links
independent observations.

There are 206,200 active models with 74,521,115,802 bytes (69.4 GiB) of
coordinate files, measured by filesystem size. All 212,549 models associated
with candidate events, including same-model cases, remain in the exported
model catalog. Original mean C-alpha confidence, fraction below pLDDT 50,
sequence hashes, versions and provenance are retained; no confidence cutoff
was silently applied. Size measurement is not raw-coordinate validation.

`metadata/duplication_model_pair_queue_completed_20260926.json` archives the
producer receipt and independent export readback. Every original event-link
column exactly matches the reviewed source tables. The readback independently
verified the complete eligible pair set, model membership, stable pair keys
and retained identical-model rows. Raw CIF byte/content validation, confidence
qualification, structural alignments, direct distance readback, domain
controls, matched nonduplication contrasts and asymmetry references remain
pending. No duplication effect is inferred from this queue.

Reproduce with:

```bash
python scripts/prepare_duplication_structure_pairs.py \
  --review results/orthology/duplication-candidate-tree-review-20260926-v1 \
  --review-plan metadata/duplication_candidate_tree_review_plan_20260926.json \
  --bridge results/structures/whole-proteome-family-coverage-20260922-v1/structure_family_bridge.sqlite \
  --catalog results/structures/whole-proteome-afdb-catalog-20260922-v1 \
  --output <fresh-output-directory>
```

Current output is
`results/structural_comparisons/duplication-model-pair-queue-20260926-v1/`.
The full protein-pair ledger and candidate guide sensitivities remain
available for downstream family- and phylogeny-aware comparisons.


## Full coordinate validation launched

The full queue model set contains 212,549 models and 76,394,885,497 bytes
(71.1 GiB) of raw CIF coordinates, including identical-model event cases.
The maximum sequence length in this frozen candidate set is 1,280 residues.
Coverage and length limits must remain explicit in subsequent comparisons.

`scripts/validate_duplication_coordinates.py` reuses the existing strict
`extract_domain_coordinates.load_atoms` reader. It verifies source hashes,
full canonical polymer sequence, atom/residue identities, duplicate atoms,
single-chain/model/alternate-location constraints, complete C-alpha coverage,
and finite coordinates, occupancy and confidence. It independently recomputes
mean C-alpha confidence and the fraction below 50 against catalog summaries.
Validated compressed records retain the full sequence, C-alpha coordinates
and per-residue confidence, plus residue counts at pLDDT >=70 and >=90.
Those counts do not yet impose an alignment eligibility threshold.

Content rejections are recorded with reasons. Missing/changed inputs detected
before validation stop a shard. Completed shards are bound to exact model
metadata and output hashes, and raw source hashes are rechecked when resuming.
An unreceipted output/partial shard requires explicit review rather than
automatic replacement. All 212,549 models must have an explicit disposition
for stage completion; completion alone is not a claim that all models passed.

Code checks use an isolated copy of a source CIF, without altering production
data. Source sequence/coordinate/confidence validation and output dimensions
passed; checkpoint reuse passed; altered confidence metadata produced an
explicit content rejection; changed raw bytes were rejected for both new
and reused checkpoints. Reproduce with:

```bash
python scripts/check_duplication_coordinate_validation.py \
  --models results/structural_comparisons/duplication-model-pair-queue-20260926-v1/models.jsonl
```

The full stage uses `metadata/duplication_coordinate_validation_plan_20260926.json`
and service `fungal-duplication-coordinate-validation-20260926.service`.
Launch identity is recorded in
`metadata/duplication_coordinate_validation_launch_20260926.json`. Four CPU
workers run with a 16-GiB memory cap and no swap. Checkpoints contain 1,000
models each (213 shards). Planning allows 32 GiB output, a 100-GiB free-disk
reserve, and an uncalibrated 0.5–24 hours. Output is
`results/structural_comparisons/duplication-coordinate-validation-20260926-v1/`.
The stage is running. Independent coordinate-output readback, confidence/PAE
sensitivity, direct alignments and biological duplication tests remain
pending. GPUs remain paused.


## Independent coordinate readback queued

The coordinate producer completed its first eight shards (8,000 accepted
models, no content rejections in those shards) before the full readback was
queued. This partial result does not establish full model-set completion.

`scripts/readback_duplication_coordinates.py` waits for the producer's exact
PID, creation time and command, then requires the completed coordinate
receipt bound to the pinned producer plan. It checks the full model-to-shard
grid and reconstructs every exported accepted C-alpha residue, sequence,
coordinate and confidence value from the raw CIF atom rows without importing
the producer extraction function. It verifies confidence summary and threshold
counts, raw/source identities, output hashes and complete disposition totals.
It shares the CIF lexical parser and residue-name dictionary with the producer;
this independence does not extend to those dependencies. Rejection identities
and reasons are checked, but rejection causes are not independently adjudicated.

The fixture matched an intact model and rejected changed coordinates,
per-residue confidence, sequence, threshold count, model identity, summary
confidence, and raw bytes. Reproduce with:

```bash
python scripts/check_duplication_coordinate_readback.py \
  --models results/structural_comparisons/duplication-model-pair-queue-20260926-v1/models.jsonl
```

The audit is queued under
`fungal-duplication-coordinate-readback-20260926.service`, using
`metadata/duplication_coordinate_readback_plan_20260926.json`. Launch identity
is in `metadata/duplication_coordinate_readback_launch_20260926.json`. Once
the producer finishes, the audit uses two CPU workers, 8 GiB RAM, no swap and
small per-shard proofs (0.1 GiB output allowance). The uncalibrated planning
range is 0.5–24 hours after the dependency completes. Output will be
`results/structural_comparisons/duplication-coordinate-readback-20260926-v1/`.
No source job is automatically restarted, and no GPU work is scheduled.
Passing readback remains a dependency for downstream coordinate use; full
structural comparisons and biological tests remain unfinished.


## Full and confidence-masked alignment inputs queued

The next stage is queued behind the exact independent-readback process. It
requires a passing audit receipt and matching coordinate producer/plan,
per-shard audit proofs and queue provenance before writing alignment inputs.
All 212,549 source models receive both `full` and `plddt70` dispositions.
Source-rejected models remain recorded. Masks with fewer than three residues
are marked too short for alignment; the native US-align source rejects
lengths below three. No failed or short input is silently replaced.

`scripts/duplication_alignment_inputs.py` serializes C-alpha-only PDBs and
retains original one-based residue numbers and the exact mask-position map.
Noncontiguous high-confidence residues keep their original numbering.
Coordinates round to 0.001 angstrom and confidence to 0.01, with field-width
and numerical-error checks. Written files are hash-verified and the manifest
records each hash, sequence, retained length and source identity.

Four native US-align fixture runs recovered the expected near-identical
scores for a known rigid rotation/translation in both input orders, using
full coordinates and a noncontiguous pLDDT>=70 mask. An independent PDB
parser checked original residue numbering; empty masks and nonfinite inputs
were tested. A synthetic completed-handoff fixture preserved six expected
full/masked/short/rejected dispositions and rejected a changed audit proof.
Run `python scripts/check_duplication_alignment_inputs.py` and
`python scripts/check_duplication_materialization.py` to reproduce.

Confidence masking is a sensitivity analysis. Masked alignment lengths and
TM-score normalizations refer to retained residues, so those scores are not
interchangeable with full-model scores. Noncontiguous masks, alignment
coverage, absent PAE qualification and interdomain-orientation uncertainty
remain explicit limitations. Matched nonduplication controls and ancestral
or outgroup references are still needed for duplication effects/asymmetry.

Plan: `metadata/duplication_alignment_input_plan_20260926.json`.
Service: `fungal-duplication-alignment-inputs-20260926.service`.
Launch identity: `metadata/duplication_alignment_input_launch_20260926.json`.
Output: `results/structural_comparisons/duplication-alignment-inputs-20260926-v1/`.
The waiting service uses one CPU, 8 GiB RAM and no swap, with a 32-GiB output
allowance and 100-GiB free-disk reserve. The uncalibrated planning interval
is 0.1–8 hours after the audit completes. At most 425,098 model/mask records
will be prepared. Production materialization and pair alignments have not
started; the coordinate producer had validated its first 20,000 models when
this dependent stage was queued.


## Full pair-alignment workload queued

`scripts/run_duplication_alignments.py` is queued behind the exact input
materializer identity. It requires the completed input receipt bound to its
plan and queue, the exact full model/mask universe, and canonical unique pair
keys. It uses the pinned US-align executable and the established monomer
options (`-mol prot -mm 0 -outfmt 0 -ter 2`).

Every one of 103,200 distinct model pairs receives both input orders under
`full` and `plddt70`: 412,800 planned directed dispositions. Native calls are
skipped only when a source input is explicitly unavailable/too short. Each
job checks input bytes, saves the complete native output and command, parses
lengths/scores/alignment sequences against the expected inputs, and saves
wall time. Exclusions, native errors, parser errors and 600-second timeouts
remain separate terminal records without silent retries or substitution.
These records must be reviewed before any inference.

Checkpoint reuse verifies plan, input-manifest, model/version, mask/order and
command bindings, source hashes and parsed metrics for completed alignments.
A bounded queue holds at most 64 pending futures. Checkpoint files retain
results, while a final manifest binds every disposition to its file hash.
Completion means every planned disposition is accounted for, not that all
alignments succeeded or passed independent geometric readback.

The runner passed a native identity-alignment fixture, cached reuse without a
native rerun, short-input exclusion, simulated native/parse/timeout failures,
and rejection of altered metrics and source coordinates. An end-to-end
completed-input fixture aligned both full-model orders and retained both
short-mask exclusions, checking all four checkpoint hashes. Reproduce with
`python scripts/check_duplication_alignment_runner.py` and
`python scripts/check_duplication_alignment_handoff.py`.

Plan: `metadata/duplication_alignment_plan_20260926.json`.
Launch: `metadata/duplication_alignment_launch_20260926.json`.
Service: `fungal-duplication-alignments-20260926.service`.
Output: `results/structural_comparisons/duplication-alignments-20260926-v1/`.
Four CPU workers run with a 16-GiB memory cap and no swap; planning allows
64 GiB output and requires a 100-GiB free-disk reserve. The uncalibrated
12–720-hour interval starts after inputs are ready and is not a guaranteed
completion bound. Per-job timings will support a measured estimate. No GPU
work or paid infrastructure is requested.

The service is waiting; production pair alignments have not begun. Native
scores and input-order differences remain descriptive pending independent
numeric readback, confidence/coverage assessment, interdomain-orientation
controls, matched ortholog/nonduplication comparisons, family/phylogenetic
modeling and ancestral/outgroup references for asymmetry.
