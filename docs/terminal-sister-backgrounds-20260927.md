# Terminal sister-pair background inventory

The duplication analysis needs comparison groups before testing an association
between duplication and structural divergence. This inventory begins that work
across all 70,307 profile-guide and 70,412 MAFFT-guide resolved gene trees.
The full inventory and independent reconstruction are complete. Guide comparison,
orthology membership and matching remain downstream; no completed matched
background or duplication-effect test is claimed.

For every internal node with exactly two terminal children, retain both genes,
their taxa, family, node and parent identifiers, reported duplication flags for
the node and parent, both terminal sequence branch lengths and their sum.
Join each gene to the frozen structure bridge and retain four availability
states: two distinct models, identical model, one model, or no models. Missing
models do not remove pairs from the inventory. Internal or multifurcating
sister groups are outside this terminal-pair design; this is not an inventory
of all possible ortholog pairs or all evolutionary events.

Each pair receives one descriptive label, with reported duplication taking
precedence: reported duplication, same-taxon unreported, or cross-taxon
unreported candidate. **An unreported duplication is not proof of speciation.**
Cross-taxon candidates can retain deeper duplication ancestry and losses.
Parent duplication flags are retained separately rather than silently filtering
them away. The two guide outputs are dependent sensitivity alternatives.

Next steps are full independent reconstruction of the inventory, cross-guide
gene-pair comparison, and ortholog-membership checks before matching candidate
backgrounds to duplicated pairs by family and sequence divergence. Domain
architecture, length, confidence, taxon sampling and shared ancestry also need
to enter eligibility, matching or downstream models. Cross-species comparisons
and same-species duplicates differ in evolutionary context; matching alone will
not justify a causal duplication interpretation. No structural outcome is used
to select this initial pool.

Run `scripts/inventory_terminal_sister_backgrounds.py --plan
metadata/terminal_sister_background_plan_20260927.json`. Inputs and imported
helpers are checksum-pinned. Output is under
`results/orthology/terminal-sister-background-inventory-20260927-v1`; each guide
has a TSV and the final receipt records exact counts and hashes. The output
directory must be new. An interrupted run must be diagnosed before any new
version is launched; no partial TSV is a completed inventory.

Known-tree fixtures in `scripts/check_terminal_sister_background_cases.py`
passed canonical gene orientation, path lengths, reported node/parent flags,
all model availability states, same-taxon and multifurcation distinctions,
multiple cherries, and rejection of invalid branches, duplicated labels and
unknown taxa. These fixtures do not replace the full production readback.

The launch reserves one CPU, 8 GiB RAM, no swap, and a 2 GiB output allowance;
the planning range is 0.1–4 hours, not a calibrated ETA. Free disk was about
12 TiB and available RAM about 942 GiB before launch. The service enforces CPU
and memory limits. No GPU predictions or paid resources are used. Process
identity and the source-bound plan are recorded in
`metadata/terminal_sister_background_launch_20260927.json`.

## Full independent readback queued

`readback_terminal_sister_backgrounds.py` waits for the recorded producer
PID/creation time/command to exit, then requires its complete checksum-bound
receipt. It traverses every source tree in postorder, independently reconstructs
parent relationships and terminal pair sets, and compares every output field
against native trees, duplication tables and the frozen model bridge. Output
order within a family is not assumed. Missing, extra or repeated pairs and
changed fields fail the check. Model coverage and candidate-status totals must
also match the producer receipt. Both readers use Bio.Phylo's Newick parser;
the checker does not invoke the producer's tree-index or pair-generation code.

The checker has its own one-CPU/8-GiB/no-swap service and a 0.1–4 hour planning
range excluding the producer wait. Its pinned config is
`metadata/terminal_sister_background_readback_plan_20260927.json`; the final
output, only on success, is
`metadata/terminal_sister_background_completed_readback_20260927.json`.
Known-tree reconstruction and deliberate corruption of five exported fields
passed in `scripts/check_terminal_sister_readback_cases.py`. These are fixture
results; the full production readback is still pending. The checker rejects
changed source or producer hashes and rechecks them after traversing all rows.

## Inventory output and guide-comparison stage

The inventory completed with 1,577,204 profile-guide and 1,577,169 MAFFT-guide
terminal pairs. Cross-taxon unreported pairs with two distinct frozen models
number 71,858 and 71,822, respectively; a further 6,446 and 6,448 map to identical
models. These counts passed full independent reconstruction; they are not counts of
eligible matched controls. The inventory includes 70,307/70,412
source trees and retains all missing-model categories.

The queued `scripts/compare_terminal_sister_guides.py` requires that full
readback to pass before proceeding. It joins the complete pair universe by
lexically ordered gene identifiers, retaining changed family/node labels,
one-guide pairs, candidate-class disagreements, parent duplication flags,
sequence distances and frozen model identities. Shared model/taxon assignments
must be identical across guides. An indexed SQLite file preserves every source
row, `guide_comparison.tsv` contains the full union, and
`modeled_candidate_union.tsv` contains pairs labeled cross-taxon unreported in
at least one guide with two available models. Identical-model pairs are retained
and distinguished from distinct-model pairs. Selection uses no structural
outcome. Shared candidate labels and unreported parents are separate flags,
not proof of orthology or independent speciation events.

The comparison plan is
`metadata/terminal_sister_guide_comparison_plan_20260927.json`; outputs will be
in `results/orthology/terminal-sister-guide-comparison-20260927-v1`. One CPU,
4 GiB RAM, no swap, and an 8 GiB disk allowance cover the disk-backed join;
0.1–4 hours is an uncalibrated planning range excluding its dependency wait.
Known-guide fixtures passed label changes, candidate disagreements, one-guide
pairs, missing/identical models, and rejection of inconsistent identities.
Production comparison and independent output checking remain pending.

## Native orthology membership queued

The next source-bound stage queries both audited native ortholog-pair streams
for every row of the modeled candidate union, retaining all positive, negative
and disagreeing assignments. Protein IDs are exact zero-based line ordinals
from the pinned `SequenceIDs.txt`, qualified with the species mapping. The ID
files are identical across guides. Candidate pairs are sorted by their native
integer keys and merged against each complete reciprocal stream. Every stream
record is checked for ordering, canonical IDs and exactly one record in each
direction; stream hashes must match the earlier full multiplicity audits before
and after querying. A negative result means absent from that native output,
not independently established paralogy. These are reconciliation-derived
ortholog assignments and are not independent biological validation.

Run `scripts/query_terminal_background_orthology.py --plan
metadata/terminal_background_orthology_plan_20260927.json`. It waits for the
recorded guide-comparison process and its complete receipt, compiles
`scripts/query_ortholog_pair_membership.cpp`, runs known-set and corruption
fixtures, and exports all candidate fields plus both membership flags under
`results/orthology/terminal-background-orthology-20260927-v1`. It can compute
before the guide-comparison independent readback, but its receipt explicitly
requires that readback and independent membership verification before any
matched-control eligibility or scientific use.

The scanner passed exact comparison to a Python set, presence and leading,
interior and trailing absence, empty streams and queries, and rejection of six
corrupted stream/query cases. The queued production scan is not yet complete.
One CPU, 4 GiB RAM, no swap and 1 GiB output are budgeted; each roughly 15.7 GiB
input stream receives a complete scan and before/after checksums. The 0.1–4 hour
planning range excludes dependency wait. No predictions or paid resources run.

## Complete guide-union readback queued

`scripts/readback_terminal_sister_guides.py` waits for the exact guide-comparison
process to finish and requires complete, mutually bound producer and inventory
readback receipts. It checks every original inventory row against its stored
SQLite JSON record and verifies table counts and database integrity. It then
merges two sorted source cursors independently of the producer's SQL union,
reconstructing every exported field and the exact modeled-candidate subset.
Missing, extra, repeated or altered rows fail exact comparison. Presence counts,
classification cross-tabs and candidate counts are recomputed from all rows.
Source and output hashes are checked before and after the full pass.

The independent merge fixtures passed shared, left-only, right-only and empty
cases, changed family/class labels and identical-model candidate flags. The
production readback is queued, not complete. Its pinned config is
`metadata/terminal_sister_guide_readback_plan_20260927.json`; successful completion
will produce `metadata/terminal_sister_guide_comparison_completed_readback_20260927.json`.
One CPU, 4 GiB RAM and no swap are enforced, with a 0.1–4 hour planning range
excluding the upstream wait. This verifies joins and classifications; it does
not establish matched comparability or a biological duplication effect.

## Full inventory verified; membership readback queued

The independent inventory reconstruction completed successfully for all
140,719 source trees and all 3,154,373 guide-specific terminal pair records.
Every pair identity, exported field, model join, sequence branch length and
summary count passed. Evidence is recorded in
[the full readback receipt](../metadata/terminal_sister_background_completed_readback_20260927.json).
The inventory and guide alternatives remain conditional on the reconciled trees.

The independent membership checker `scripts/readback_background_orthology.py`
waits for both the membership producer and full guide-union checker. It requires
both complete, bound receipts, independently rebuilds protein ordinals and all
copied candidate fields, then uses fixed-record binary searches for every
candidate in each native stream. This differs from the producer's sequential
merge. Search records must contain reciprocal entries; absent queries require
appropriate bracketing. The prior full stream audit and before/after hashes
establish global ordering and multiplicity outside the accessed records.

Known-set checks passed all 91 possible pairs in a 14-ID universe, plus empty
stream and corrupt-record cases. The queued production check uses one CPU,
4 GiB RAM, no swap and a 0.1–4 hour planning range excluding its waits. Config:
`metadata/terminal_background_orthology_readback_plan_20260927.json`; success
output: `metadata/terminal_background_orthology_completed_readback_20260927.json`.
This validates native output membership, not independent biological orthology,
matched comparability or a causal effect of duplication.

## Guide comparison produced; measurement inventory queued

The guide comparison produced a union of 1,582,382 terminal pairs: 1,571,991
shared, 5,213 profile-only and 5,178 MAFFT-only. The modeled candidate union
contains 78,372 pairs, including 78,202 shared candidates, 102 profile-only and
68 MAFFT-only. Of the shared modeled candidates, 70,684 have unreported parents
in both guides. Shared pairs do not change the reported-duplication versus
cross-taxon-unreported label between guides. These are producer results;
independent guide-union readback and orthology membership remain running.

`prepare_background_measurements.py` is queued behind successful full native
membership readback, which itself requires the guide-union readback. It retains
all 78,372 candidate rows and adds separate per-guide flags requiring both the
cross-taxon-unreported label and native ortholog membership. Either-guide,
both-guide and both-guide/unreported-parent flags remain separate sensitivity
sets. These flags are computational qualification, not matched comparability.
Pairs qualifying in neither guide remain explicit in the output ledger.

Every gene is joined to the frozen sequence-qualified model bridge and catalog.
Model IDs, versions and coordinate paths must agree. The ledger adds sequence
hashes, lengths, mean C-alpha pLDDT and the fraction below 50 for each model.
Distinct model pairs qualified in either guide are partitioned into already
present in the primary queue, already present in the reference queue, or new.
Existing queue membership does not establish completed or valid measurements.
Identical models remain explicit without redundant alignment jobs. All candidate
models, active models and additional models absent from both prior model
inventories are exported separately; additional does not mean newly predicted.

Plan: `metadata/background_measurement_inventory_plan_20260927.json`.
Output: `results/structural_comparisons/background-measurement-inventory-20260927-v1`.
Known-case eligibility, pair-symmetry and version-sensitivity checks passed.
One CPU, 4 GiB RAM, no swap, 2 GiB output and a 0.1–2 hour planning range excluding
waits are budgeted. This stage parses no coordinates and launches no predictions
or alignments. Full output readback, additional-coordinate validation,
domain/coverage comparisons and matching remain necessary.

## Guide union verified and pre-matching support queued

The full independent guide-union check passed all 3,154,373 source records,
1,582,382 union rows and the exact 78,372-row modeled subset. See
[readback evidence](../metadata/terminal_sister_guide_comparison_completed_readback_20260927.json).
Native membership production also completed: 77,863 present/509 absent in the
profile stream and 77,857 present/515 absent in MAFFT. Independent membership
validation is still running; absence is not itself proof of paralogy.

`assess_background_matching_support.py` is queued after successful membership
readback. It retains all 218,473 reviewed guide-specific two-model duplicate
targets, including identical-model targets. The target pair, family, taxon and
node must match the independently audited terminal inventory. For each guide,
three background sets are considered: that guide's cross-taxon candidates with
native orthology, candidates qualified in both guides, and the latter with
unreported parents in both guides. All available models, including identical
models, remain in these feasibility pools.

For each target and set, count same-family backgrounds, those involving the
focal taxon, and those involving neither focal taxon. Also count each category
within multiplicative sequence-distance ranges of 1.25, 1.5 and 2.0 around the
target. These are descriptive sensitivity ranges, not calibrated matching
criteria. A zero target distance matches only exactly zero background distances;
no pseudocount is introduced. Pool counts are not numbers of independent events
or guarantees of usable structural comparisons. Targets with no support remain
in the 655,419-row output. No outcome-dependent target selection is performed.

Plan: `metadata/background_matching_support_plan_20260927.json`; output:
`results/orthology/background-matching-support-20260927-v1`. One CPU, 4 GiB RAM,
no swap, 1 GiB output and 0.1–2 hours excluding wait are budgeted. Per-family
sorted arrays and binary interval counts avoid Cartesian pair expansion.
Eighteen enumerated interval fixtures, zero/empty cases and nested qualification
sets passed. Production support counts and full independent validation remain
pending. This diagnoses possible support only: domain architecture, length,
confidence, coverage, taxonomic balance and shared ancestry still require
assessment before any matched effect test.
