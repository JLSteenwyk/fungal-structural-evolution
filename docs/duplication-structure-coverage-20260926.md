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
