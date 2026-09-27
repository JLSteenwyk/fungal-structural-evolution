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
