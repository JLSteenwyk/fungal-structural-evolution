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

## Membership and measurement inventories verified; coordinates running

All 78,372 membership rows passed independent binary-search verification,
including exact protein ordinals and absent records. The guide comparison and
membership receipts are now both complete. The measurement inventory also passed
full independent reconstruction of every field, frozen model join, metadata
value, pair identity and model/work partition. Evidence:
[orthology readback](../metadata/terminal_background_orthology_completed_readback_20260927.json)
and [measurement readback](../metadata/background_measurement_inventory_completed_readback_20260927.json).

Of 78,372 candidates, 77,909 meet candidate-plus-native-ortholog criteria in at
least one guide, 77,707 in both, and 70,209 also have unreported parents in both.
The ledger retains 463 candidates qualified in neither guide and 6,448 identical-
model pairs. There are 71,461 qualified distinct-model pairs: 11 already in the
reference queue and 71,450 new pairs. The candidate universe uses 150,280 models,
of which 149,356 are active for the qualified pool and 148,104 are absent from
the prior primary/reference model inventories. These are existing predictions,
not a requirement to infer 148,104 new structures.

`validate_background_coordinates.py` is now validating every additional model,
using the existing pinned raw-CIF/shard validator. The input contains 66,498,349
residues and 63,311,748,521 coordinate bytes. Four CPU workers, 16 GiB RAM, no
swap, a 32 GiB output allowance and a 100 GiB free-disk reserve are configured.
The 0.5–24 hour interval is an uncalibrated planning range. Checkpoints contain
1,000 models and require unchanged source hashes for reuse. Content rejections
remain explicit; missing or changed files fail provenance checks. No new
predictions or structural alignments are launched by this stage.

Plan: `metadata/background_coordinate_validation_plan_20260927.json`; output:
`results/structural_comparisons/background-coordinate-validation-20260927-v1`.
The coordinate-array/checkpoint/rejection fixture passed with the additional
model inventory. An initial fixture invocation omitted its required `--models`
argument and exited with usage output; the corrected invocation passed before
launch. Full independent raw-coordinate readback remains downstream, followed
by input materialization, background alignments and domain/coverage matching.

The matching-support producer also completed all 655,419 target/set rows, but
that independent readback remains pending. Its provisional profile-guide counts
show 31,435 targets with no same-family background under the local-guide set;
only 53,184 have any background within the factor-1.5 distance range. These
preliminary support counts motivate retaining unsupported targets and reporting
the eventual matched estimand separately from the full target universe.

## Matching-support counts independently verified

Full independent readback passed all 218,473 target records and all 655,419
qualification-set rows. The checker reconstructs original target identities and
qualified background pools, then uses unsorted vector masks rather than the
producer's sorted binary interval counts. Every family, focal/nonfocal and
sequence-window count, the complete target/set grid and all summaries matched.
Evidence: [full support readback](../metadata/background_matching_support_completed_readback_20260927.json).

| Target/support category | Profile guide | MAFFT guide |
| --- | ---: | ---: |
| All reviewed duplicate targets | 109,245 | 109,228 |
| No same-family background, local-guide qualification | 31,435 | 31,414 |
| At least one background within factor 1.5, local-guide qualification | 53,184 | 53,188 |
| Within factor 1.5, both guides and unreported parents | 48,571 | 48,566 |
| Same as preceding row, with focal taxon represented in background pair | 5,103 | 5,100 |

Thus even sequence-distance/family support covers only part of the original
target universe, and focal-taxon representation is much more limited. These
are available-pool counts before length, domain, confidence and alignment
filters, not final match counts. Unsupported targets must remain in the
coverage report; an effect estimated in a selected matched subset cannot be
silently generalized to every duplicated gene or fungal lineage. Both guide
columns describe dependent sensitivity analyses, not independent replication.
The factor-1.25 and factor-2 alternatives and exact-zero cases remain in the
complete table and receipt.

## Complete background-coordinate readback queued

`readback_background_coordinates.py` waits for the exact coordinate producer,
requires its complete receipt and exact model/shard mapping, then uses the
pinned independent raw-CIF checker for every exported record. Accepted sequence,
C-alpha positions, pLDDT values and summary counts are reconstructed from atom
rows. Rejection identities and reasons are retained, but rejection causes are
not independently adjudicated by this checker. All raw/shard hashes are checked.
The wrapper adapts the existing full-coordinate checker to the background
producer status without changing any scripts pinned by existing jobs.

Plan: `metadata/background_coordinate_readback_plan_20260927.json`; output:
`results/structural_comparisons/background-coordinate-readback-20260927-v1`.
Four CPUs, 16 GiB RAM, no swap and 0.1 GiB output are budgeted, with an
uncalibrated 0.5–24 hour interval excluding producer wait. The fixture accepted
exact raw reconstruction and rejected six altered exports plus changed raw
bytes. Production coordinate validation is running; its full readback is queued.

## Background domain controls

The background models and pairs now have the same four-policy domain annotation
inventory as the duplicate targets. `inventory_background_domain_controls.py`
requires the independently verified measurement inventory and frozen Pfam
registry. It joins all 150,280 candidate models by exact sequence, model version,
length and path, then compares every one of the 71,461 qualified distinct-model
pairs under alignment/envelope and E-value/bitscore policy alternatives.

Ordered annotation signatures retain all Pfam types and repeated occurrences.
The five categories distinguish same ordered annotations, same content in a
different order, different content, one unannotated model, and neither annotated.
Single-copy Domain matches require one occurrence of that versioned accession
in each model and conservative interval eligibility on both sides. A second
ineligible copy still blocks single-copy classification. Missing annotations
are not biological absence, and annotation differences are not inferred gains,
losses or rearrangements. Residue-level confidence/coverage and domain-coordinate
comparison remain downstream.

Producer output contains all model annotations, 285,844 pair/policy records and
every matched single-copy domain interval pair. Under alignment/E-value policy,
50,539 pairs have the same ordered annotations, 5,867 differ in content, 11
share content in a different order, 2,833 have only one annotated model and
12,211 have neither annotated. There are 36,717 matched domain pairs under that
policy; the other three policies yield 36,683, 36,580 and 36,549. All counts passed full independent readback and are not independent
evolutionary events or matched-target counts. See the
[complete domain readback](../metadata/background_domain_control_completed_readback_20260927.json).

Plans are `metadata/background_domain_control_plan_20260927.json` and
`metadata/background_domain_control_readback_plan_20260927.json`; data are under
`results/structural_comparisons/background-domain-controls-20260927-v1`.
Both stages use one CPU, 16 GiB RAM and no swap; output allowances are 2 GiB and
0.01 GiB, respectively, with uncalibrated 0.1–6 hour planning ranges. No new
predictions or native structural alignments are launched. The readback joins
segments and policy membership separately and independently reconstructs every
exported annotation, architecture record and domain-pair match. Known cases and
1,000 independent pair fixtures passed.

The producer completed before its live PID could be captured. That prevented
creation of the first readback plan, and the initial readback service exited
with a missing-plan error before reading data. After verifying the completed
producer receipt and artifact hashes, a receipt-bound readback plan was created
and a new `fungal-background-domain-readback-20260927-v2` service launched.
The failed launch and replacement are recorded in
`metadata/background_domain_control_produced_20260927.json`; no producer data
were rerun or overwritten.

## Conserved-architecture support assessment running

`assess_architecture_matched_support.py` extends the verified family/sequence
support diagnostic to all four domain policies. Every one of the 218,473 target
records is retained under all three native-orthology qualification sets and all
four policies, producing 2,621,676 target/set/policy records. A qualified pool
requires the two background proteins and both target proteins to have the same
nonempty ordered signature of versioned Pfam accessions and types, with every
retained hit marked conservative under that policy. Repeats, order and versions
remain part of identity. This defines a within-conserved-architecture analysis;
it does not replace the separate aim of studying architecture-changing events.

Targets with neither model annotated, one model unannotated, different ordered
annotations, or a shared but nonconservative signature retain explicit status
and zero support under this particular comparison rule. Missing annotations are
not coded as a matching empty architecture. Identical-model targets/backgrounds
remain represented. Focal/nonfocal counts and all three sequence-distance ranges
are recomputed within each family/signature pool and must not exceed the prior
unrestricted support counts. The analysis uses no structural response variable.

Plan: `metadata/architecture_matching_support_plan_20260927.json`; output:
`results/orthology/architecture-matched-background-support-20260927-v1`.
One CPU, 8 GiB RAM, no swap and 4 GiB output are budgeted with a 0.1–4 hour
planning range. Known shared/missing/nonconservative/content/repeat/order cases
passed (`scripts/check_architecture_support_cases.py`). Production is running;
full independent output validation remains required. Length, prediction
confidence, coordinate coverage and phylogenetic dependence are still outside
this support diagnostic and remain necessary for final matching/inference.

## Full architecture-support readback running

The producer completed all 2,621,676 records. An independent readback now rebuilds
all target identities, signatures, qualification pools and architecture labels,
then recounts each family/focal/nonfocal distance window with unsorted NumPy
masks rather than binary interval searches. The checker verifies the exact full
row grid and every summary count; it samples no rows. One thousand independent
signature/status fixtures and an empty-annotation guard passed.

Config: `metadata/architecture_matching_support_readback_plan_20260927.json`;
success output: `metadata/architecture_matching_support_completed_readback_20260927.json`.
One CPU, 8 GiB RAM, no swap and 0.01 GiB output are budgeted, with an uncalibrated
0.1–4 hour planning range. Production verification remains running.

Provisional producer counts under alignment/E-value policy and the strict
both-guide/unreported-parent qualification set show 18,454 targets per guide
with at least one same-family, same-conservative-architecture background within
factor 1.5 sequence distance; 1,322 also have focal-taxon representation. The
profile guide has 34,819 shared-conservative-architecture targets overall,
13,100 shared but nonconservative, 10,192 different ordered annotations, 8,591
with one model unannotated and 42,543 with neither annotated. The full readback
must pass before treating these counts as verified. These distinctions preserve
architecture-changing and annotation-uncertain targets for separate analysis;
this conserved-architecture subset cannot stand for the entire project.

## Architecture support verified; alignment inputs queued

The full independent architecture-support readback passed all 218,473 targets
and 2,621,676 target/set/policy records. The previously provisional counts above
are now verified, including 18,454 targets per guide with strict-set,
alignment/E-value, factor-1.5 architecture-compatible support and 1,322 with
focal-taxon representation. Evidence:
[complete architecture-support readback](../metadata/architecture_matching_support_completed_readback_20260927.json).
This does not establish final matches or extend inference to unsupported targets.

`materialize_background_inputs.py` is queued behind the exact background
coordinate-readback process. It requires a complete bound coordinate audit,
validated model/shard identities, unchanged source and proof hashes, and the
verified additional-model partition. All 148,104 additional models receive
full and pLDDT≥70 dispositions (296,208 expected records). Original sequence
positions are retained; rejected content and masks with fewer than three
residues remain explicit. Existing primary/reference inputs for the other
background models will be joined separately with their existing proofs.

Serialization uses the existing pinned C-alpha writer. PDB coordinates are
rounded to 0.001 Å, confidence to 0.01, with overflow and rounding checks.
Masked comparisons normalize to retained length and must not be conflated with
full-model scores. This stage does not apply PAE screening, infer structures,
or run alignments. A fresh output directory prevents overwriting prior work.

Plan: `metadata/background_alignment_input_plan_20260927.json`; output:
`results/structural_comparisons/background-alignment-inputs-20260927-v1`.
One CPU, 8 GiB RAM, no swap, 32 GiB output and a 100 GiB free-disk reserve are
specified; 0.5–24 hours is an uncalibrated planning range excluding the wait.
The residue count gives an upper estimate near 10.8 GB of PDB atom text before
metadata and filesystem overhead. The complete synthetic handoff passed all
six fixture dispositions, sparse residue mapping, written hashes and changed
proof rejection. Four native rigid-transform alignments and full/noncontiguous
mask serialization checks also passed. Full production input readback and
background comparisons remain downstream.

## Background structural comparisons queued

Queued `fungal-background-alignments-20260927.service` with versioned plan and
launch identity. It waits for background input materialization, then checks the
complete, disjoint primary/reference/background input partitions, coordinate
readback bindings, exact raw-model provenance, both-mask disposition grids and
background inventory proof. All 71,450 new model pairs receive full and pLDDT70
comparisons in both input orders (285,800 maximum native calls/dispositions).
The 11 existing reference pairs are excluded from computation and must later be
joined from independently checked reference results. Identical-model and
ineligible candidates remain explicitly recorded in the original inventory.

Two CPU workers, 8 GiB RAM, no swap, 600-second per-call timeout and a 100 GiB
free-disk reserve are enforced. The output estimate is 48 GiB; the uncalibrated
12–1,050-hour planning interval scales the reference workload and excludes
waiting. This is not an ETA or a guaranteed bound. No GPU prediction or paid
resources were enabled. Ready input bytes are hash-checked before and after
each native comparison; excluded inputs, errors and timeouts remain explicit,
with no automatic substitution or retry.

The three-source fixture passed four native alignments and four short-mask
exclusions, existing-pair reuse, all checkpoint hashes, and rejection of altered
coordinate proofs, raw provenance and missing masks. Production measurements
remain pending. Independent numerical readback must precede scientific use;
matched effect estimation, domain/orientation controls and shared-ancestry
modeling remain outstanding. This stage measures whole-chain background
comparisons and does not establish biological orthology or duplication effects.

## Full background numerical readback queued

Queued an independent checker behind the background alignment controller. It
reconstructs the exact new-pair/two-mask/two-order grid from the frozen inventory,
checks every checkpoint and all source bindings, then independently parses PDB
residue positions, native alignment text, sequence identities and least-squares
RMSDs for every successful comparison. It does not invoke the producer handoff
or parser. The existing numerical helper remains unchanged because primary and
reference audits also pin it. TM-scores are checked against native text, not
independently reoptimized; errors/exclusions are checked for consistency rather
than independently rerun. All production numerical results remain pending.

The complete fixture passed eight dispositions including four native numeric
reconstructions. Deliberately altered RMSD metrics, with all enclosing file
hashes updated, were rejected. Production retains the strict 0.00501 Å printed
RMSD tolerance. The earlier domain audit's small-core discrepancy is not waived:
if a similar discrepancy appears here, this audit stops for explicit review.

Resources: one CPU, 16 GiB RAM, no swap, estimated 1 GiB output and uncalibrated
1–48 hours after producer completion. All 285,800 possible dispositions are in
scope. No GPU or paid infrastructure is used. Plan, script hashes, dependency
identity and checker launch are versioned in metadata.

## Taxonomic concentration of matching support

Regrouped all 2,621,676 verified architecture-support rows into 3,672
separate guide/background-set/policy/taxon cells. The output preserves all
153 observed target taxa, including zero-support taxa, architecture categories,
identical-model target counts, unique family counts and support at all three
sequence-distance ranges, with focal-taxon counts separately. An independent
pandas implementation verified every aggregate and concentration statistic.
Commands, hashes and all 24 guide/set/policy summaries are recorded in
`metadata/architecture_support_taxa_completed_20260927.json`; the full table is
`results/orthology/architecture-support-taxa-20260927-v1/taxon_support.tsv`.

For the alignment-E-value policy and the stricter both-guide/unreported-parent
background set, a multiplicative 1.5 sequence-distance range yields 18,454
supported target records in either guide. They occur in 143 of the 153 observed
target taxa; only 81 taxa have any supported target with a background that
includes the focal taxon. Five taxa account for 36.469% of supported targets:
F1507870, F5486, F28583, F5219 and F423460. Guide totals are sensitivity results,
not independent replicates or counts to add together.

This is availability, not completed matching. The denominator is modeled
terminal duplicate targets, not the complete 526-taxon sampling or all genes.
Taxa absent from that observed target universe are not represented by zero rows.
Identical-model targets remain explicit, and missing/different architectures
remain outside this conserved-architecture comparison without being discarded
from the larger project. Counts do not establish independent evolutionary events.

The concentration motivates reporting taxon-balanced as well as event-weighted
effects, leave-one-taxon-out sensitivity, and family/shared-ancestry adjustments
when structural responses are ready. Equal taxon weighting alone does not
correct phylogenetic dependence, missing structural coverage or assembly and
gene-copy uncertainty. Domain-architecture-changing events remain a separate
analysis, rather than being generalized from these conserved controls.

## Background domain interval inventory verified

Expanded all 146,529 single-copy Domain match records across the four annotation
policies into both alignment and envelope boundaries, preserving every original
model-pair, policy, accession and hit link. The resulting 293,058 links refer to
131,986 distinct intervals and 66,929 interval pairs across 54,004 models.
Intervals contain 21,223,719 residues in total, with a maximum length of 1,046.
This defines at most 267,716 directed full/pLDDT70 comparison dispositions
before checking for reusable target-domain computations or short masks.

Independent full reconstruction passed every source field, interval identity,
bound, sequence/raw-model provenance and exact deduplicated interval/pair set.
The interval helper also rejected invalid bounds and distinguished model
versions. Reproducible commands, source hashes and completion/readback records
are in `metadata/background_domain_pair*20260927.json`; large tables remain in
`results/structural_comparisons/background-domain-pairs-20260927-v1`.
The metadata expansion and readback were planned for one CPU, 8 GiB RAM and
2 GiB output (producer), with an uncalibrated 0.02–2-hour range per stage.

No new domain alignment is claimed. Next steps are reuse partitioning against
existing domain computations, coordinate input materialization after full raw
coordinate readback, both-mask/order comparisons and independent numeric
verification. Within-domain results will help separate fold changes from
whole-chain orientation differences. They do not alone establish domain
homology, biological orthology, inherited change, or functional effects.

## Domain overlap checked and coordinate extraction queued

Exact inventory comparison found 737 shared interval descriptors but **zero
shared domain pairs** between the background inventory and the existing target/
reference domain inventory. All shared descriptors, including raw-source hashes,
sequence identities and bounds, agree. Thus all 66,929 background domain pairs
require new measurements; no earlier pair result is being reused. The overlap
record does not certify existing coordinate serialization or native results.
Reproduce it with `scripts/check_background_domain_overlap.py --background
results/structural_comparisons/background-domain-pairs-20260927-v1 --existing
results/structural_comparisons/duplication-domain-pair-inventory-20260926-v1
--output metadata/background_domain_overlap_20260927.json`.

Queued `fungal-background-domain-inputs-20260927.service` behind the verified
live background coordinate-audit identity. It streams the three disjoint audited
coordinate collections and extracts all 131,986 domain intervals across 54,004
source models, producing 263,972 full/pLDDT70 dispositions. Every shard proof,
source identity and complete interval universe is checked. Original residue
positions are retained, with rejected sources and short masks explicit.

The three-source six-disposition fixture passed independent PDB coordinate,
confidence and original-position parsing, invalid-bound checks and altered-proof
rejection. One CPU, 8 GiB RAM, no swap, 8 GiB estimated output and a 100 GiB disk
reserve are configured. Planning is 0.5–24 hours after the audit wait; this is
uncalibrated. No GPU or paid resources are used. Production extraction,
independent serialization readback and domain comparisons remain pending.

## Full domain serialization readback queued

Queued `fungal-background-domain-input-readback-20260927.service` behind the
verified live extraction controller. The checker reconstructs every interval's
full/pLDDT70 residue selection independently, then parses every ready PDB atom
for original residue number, amino acid, XYZ, confidence and occupancy. It checks
the exact 131,986-interval/two-mask grid, all source identities and shard proofs,
rejected and short-mask dispositions, counts and output bytes. Coordinate
receipt bindings are checked again after the complete readback.

The three-source fixture passed six dispositions and all 16 emitted CA atoms.
A deliberately changed coordinate failed numerical reconstruction despite an
updated PDB hash. Source helpers pinned by other jobs remain unchanged. The
checker does not repeat the raw-CIF audit or validate biological domain bounds.

Resources: one CPU, 8 GiB RAM, no swap, negligible receipt output and an
uncalibrated 0.5–24 hours after extraction. At most 42.45 million atom visits
precede the effect of confidence masks. Production serialization verification
and all background domain alignments remain pending; no GPU was enabled.

## Full background domain alignments queued

Queued `fungal-background-domain-alignments-20260927.service` behind the verified
live domain serialization checker. The handoff requires the complete independent
input proof, exact inventory, complete two-mask grid and consistent interval
identities. All 66,929 interval pairs receive both input orders and full/pLDDT70
masks, for up to 267,716 native calls/dispositions. Native output and timing,
input hashes, excluded inputs, errors and timeouts are checkpointed explicitly.
Independent numerical reconstruction remains required before scientific use.

Native handoff fixtures passed both orders, short-mask exclusions, checkpoint
reuse/hashes and altered-proof rejection. The immutable native worker used by
prior domain runs is reused. No GPU prediction is enabled.

Resources: two CPU workers, 8 GiB RAM, no swap, 24 GiB output allowance,
100 GiB free-disk reserve and 600 seconds per native call. Sampling every 1,000th
prior target-domain checkpoint gave 282 full/order-zero dispositions with mean
0.06168005 seconds and median 0.02650365 seconds. A mechanical two-worker
extrapolation is 2.29 hours but omits hashing, I/O, audit, contention and shifts
in domain composition, and does not sample other masks/orders. The broader
4–72-hour planning range excludes dependencies and is not a guaranteed bound.
All production alignments and numeric verification remain pending.

## Full background domain numerical verification queued

Queued `fungal-background-domain-alignment-readback-20260927.service` behind
the verified live domain alignment controller. It checks the complete 267,716
pair/mask/order disposition grid, checkpoint hashes, source identities and
commands; every successful residue correspondence, exact sequence identity and
least-squares RMSD is independently reconstructed from hashed PDBs. Native
TM-scores are checked against the native text, not independently optimized.
Excluded/error/timeout cases are checked for consistency, not rerun.

The native fixture passed both input orders and excluded masks; deliberately
omitting a disposition was rejected despite an updated manifest hash. The
strict 0.00501 Å RMSD printed-rounding tolerance is unchanged. A discrepancy
stops this audit with the precise pair/mask/order and checkpoint path in the
error, preserving failures for review rather than qualifying them as passed.
This matters given the previously recorded small-core discrepancy in the target
domain audit; that earlier failure is neither erased nor waived here.

Resources: one CPU, 8 GiB RAM, no swap, estimated 1 GiB output and uncalibrated
1–48 hours after alignments finish. The checker performs no native reruns and
uses no GPU or paid infrastructure. Production measurements, numerical
verification, boundary/confidence sensitivity and biological effects remain
pending. The complete background domain pipeline is now queued through its
numerical verification stage.

## Control availability figure

![Verified background control availability and taxonomic concentration](figures/background_support_coverage_20260927.png)

[Vector PDF](figures/background_support_coverage_20260927.pdf) and
[SVG](figures/background_support_coverage_20260927.svg) are available, with
[coverage points](figures/background_support_coverage_20260927_coverage.tsv) and
[concentration points](figures/background_support_coverage_20260927_concentration.tsv).

Panel A reports the percentage of all observed modeled terminal duplicate
records with at least one same-family, conservative shared-architecture
background within a factor of 1.5 of sequence distance. Backgrounds require
native ortholog assignment and unreported parent duplication in both guides;
this does not establish speciation. Filled points allow either background taxon;
open points require the focal duplicate taxon among background endpoints.
Approximately 17% have any such control, compared with about 1.2% with focal
support. All four annotation policies and both guide alternatives are shown.

Panel B ranks all 153 observed target taxa by supported-target count under the
alignment-E-value policy. Five taxa contribute 36.5% of supported targets.
The dotted diagonal is an equal-share reference, not a statistical null test.
Identical-model targets and zero-support observed taxa remain explicit; taxa
outside the observed target universe are outside these denominators. Guides are
sensitivity alternatives, not independent replicates. This figure measures
availability, not matched effects, uncertainty intervals or causal inference.

Every one of the 16 coverage and 306 concentration points was independently
reconstructed from verified source rows, including denominators, ranks and
cumulative counts. The rendered PNG was visually inspected. Reproduction
commands and script/export hashes are recorded in
`metadata/background_support_figure_reproduction_20260927.json`.

## Full candidate matching graph launched

Started `fungal-background-match-graph-20260927.service` to enumerate actual
candidate controls, extending the verified support counts. All same-family,
conservative shared-architecture target/background pairs within a factor of two
of sequence distance are included for both guides and all four annotation
policies. Exact zero distances match only zero; tighter 1.25/1.5 distance flags
and focal-taxon membership are retained. Background nodes preserve both-guide
native-orthology and unreported-parent flags, so all three qualification sets
can be reconstructed without duplicating graph edges.

The full verified support table predicts 1,661,948 edges. All 218,473 target
records and 873,892 target/policy dispositions remain, including unsupported
ones. Target and background nodes retain gene/taxon/model/version identities,
raw and sequence hashes, endpoint lengths, mean CA pLDDT and low-confidence
fractions. These covariates are preserved for balancing, not yet used to claim
that a control is well matched. Structural responses do not enter construction.

The indexed distance search passed 1,000 direct-enumeration fixtures, exact-zero
and boundary cases, empty pools and invalid distances. Production cross-checks
every emitted target/policy's distance and focal counts against the independently
verified support table. Full independent graph reconstruction remains pending.
Reused controls and shared genes, taxa, families or guide alternatives must not
be treated as independent observations. Final control selection, covariate
balance, outcome eligibility, phylogenetic adjustment and effect estimation
remain outstanding.

Resource plan: one CPU, 16 GiB RAM, no swap, 4 GiB output allowance and
uncalibrated 0.1–4 hours. This is metadata processing only; no GPU prediction,
native structural comparisons or paid infrastructure was started.

## Full candidate graph produced; independent verification running

Production completed all 1,661,948 expected edges, 218,473 target nodes,
155,616 guide-specific background nodes and 873,892 target/policy dispositions.
The producer exited successfully. Started an independent full-data checker that
rejoins every node field to source genes, model descriptors and guide distances;
it does not import graph construction or annotation-classification helpers.

Every unique edge must satisfy the exact guide, family, conservative architecture
and sequence-distance requirements, with correct focal and tighter-range flags.
The checker reconstructs counts for all three background sets and all three
distance ranges and compares every one of 2,621,676 independently audited support
rows. Valid unique edges plus equal eligible counts establish completeness of
the graph; all unsupported target/policy dispositions are checked too. Source
annotation-status cases passed before launch. Production verification remains
pending, so the graph is not yet qualified for selecting controls.

Resources: one CPU, 24 GiB RAM, no swap, negligible receipt output, uncalibrated
0.1–4-hour planning range; no native structural comparisons, GPU or new charges.
The full objective, including architecture-changing events and phylogenetic
adjustment, is unchanged by this conserved-architecture candidate graph.

## Candidate graph verified; endpoint covariates running

Full independent graph verification passed all 218,473 target nodes, 155,616
background nodes, 1,661,948 edges, 873,892 target/policy dispositions and all
2,621,676 support-table rows. The proof is archived as
`metadata/background_match_graph_completed_readback_20260927.json`.

Started full-edge length/confidence characterization. For each of the two
possible endpoint correspondences, retain the maximum endpoint length ratio,
absolute mean-CA-pLDDT difference and absolute low-confidence-fraction difference.
Three prespecified descriptive sensitivity bands require all three quantities
to be within (1.1, 5, 0.05), (1.25, 10, 0.1) or (1.5, 15, 0.2), respectively.
At least one *single* correspondence must pass all components; favorable values
from different correspondences cannot be mixed. These are matching tolerances,
not validated model-reliability thresholds. Both raw mappings remain available.

All graph edges remain explicit. Shared genes, model/version identities and
sequence hashes are counted to expose dependence and potential trivial matches.
No exclusion or final match selection is claimed. Confidence balancing can
change the estimand and cannot eliminate sequence-derived prediction circularity.
Length/confidence handling does not replace alignment coverage, PAE, domain
orientation or phylogenetic controls.

Known cases passed identity, endpoint swaps, invalid covariates and a deliberate
case where lengths match only under one orientation and confidence only under
the other (correctly failing joint calipers). Resource plan: one CPU, 8 GiB RAM,
no swap, 2 GiB output, uncalibrated 0.1–4 hours. Production and full independent
covariate verification remain pending; no structural responses or GPUs are used.

## Covariates independently verified

Completed metadata differences and dependence flags for every one of the
1,661,948 graph edges. Started independent reconstruction using vector arrays
in 50,000-edge chunks rather than the producer's scalar endpoint loop. It checks
every original graph field and row order, both sets of numeric differences,
all joint tolerance flags, distinct shared gene/model/sequence counts and every
guide/policy summary. Exact node floats are loaded through JSON parsing; TSV
floats use round-trip parsing to avoid introducing precision loss in the audit.

One thousand independent scalar/array fixtures passed, along with invariance of
eligibility under endpoint swapping. Source proof and all completed artifacts
are bound by hash. Full production readback passed all 1,661,948 edges. One CPU, 8 GiB RAM,
no swap and negligible receipt output are allocated; the uncalibrated planning
range is 0.1–4 hours. No structural outcomes enter this check and no matches or
biological effects have been selected or estimated.

The checker completed before its live PID could be captured; its completed proof
and inactive/dead service with exit status zero were verified. It was not
restarted. The proof and complete verified production summary are archived in
`metadata/background_match_covariate_completed_readback_20260927.json` and
`metadata/background_match_covariates_completed_20260927.json`.

For profile/alignment-E-value, 55,785 of 207,950 candidate edges meet the tight
band, 127,538 the moderate band and 173,094 the wide band. No edge in this subset
shares a gene, model or sequence with its target. These are dependent edge
counts, not unique supported targets, selected matches or effective sample sizes.

## Metadata control selection launched

Started deterministic nearest-control selection over the full verified graph.
The sensitivity grid has 54 scenarios: three native-background qualification
sets × three sequence-distance ranges × three metadata tolerance bands × any-
or focal-taxon backgrounds. Both guides and all four annotation policies remain
separate. Every one of 873,892 target/policy records retains explicit matched
and unmatched scenario lists: 47,190,168 scenario decisions. Selected records
are sparse; the complete target universe and all denominators remain intact.

Within an eligible scenario, the fixed ranking score sums squares of:

- log(background/target sequence distance) divided by log(1.5);
- log(maximum endpoint length ratio) divided by log(1.25);
- maximum endpoint mean-pLDDT difference divided by 10;
- maximum endpoint low-confidence-fraction difference divided by 0.1.

A single endpoint correspondence must pass all calipers; among eligible
correspondences the lowest score is used. Zero target sequence distance admits
only zero background distance, contributing zero to that score component.
Exact ties use background node hash and then endpoint order. These scales are
fixed design choices, not learned effect-size weights or reliability guarantees.
No structural response is used to construct or rank controls.

Selection is one control per target/scenario with replacement. Shared target/
background genes, model identities or sequence hashes are excluded explicitly.
Identical-model pairs within a target or within a background remain marked in
the source node records. Reuse counts are exported separately for each guide,
policy and scenario. Matching does not make repeated controls, related taxa,
gene families or alternative guides independent. Outcome eligibility and
confidence adjustment can change the estimand and require explicit accounting.

Known tests passed all scenarios, ties, shared/focal/qualification exclusions,
empty pools, exact zero, distance boundaries and the joint-orientation
counterexample. Resources: one CPU, 16 GiB RAM, no swap, 8 GiB output allowance,
uncalibrated 0.25–8-hour planning range. Production, full independent selection
readback, measured covariate balance and phylogenetically informed effect tests
remain pending. No GPU or paid resources were enabled.

## Selection produced; full independent reconstruction running

Production completed 2,786,912 selected target/control records across 54
scenarios, with all 873,892 target/policy records and their unmatched scenarios
retained. These are sensitivity-grid records, not independent pairs to pool.
Started an independent checker that enumerates eligible backgrounds and both
endpoint mappings directly for every one of the 47,190,168 scenario decisions.
It reconstructs chosen identity/order/score, eligible and tied-candidate counts,
runner-up score gaps, unmatched lists, all summaries and per-scenario reuse.
It imports no selection or candidate-ranking helpers from the producer.

The initial randomized fixture mistakenly serialized a shared-identity Boolean
as `True`/`False` instead of numeric `0`/`1`, so it failed before checking choices.
The shell continued to launch the production checker. Its pinned code/plan were
left unchanged; a separate corrected fixture passed 5,400 randomized scenario
comparisons plus all 54 exact-tie cases, including zero distances and exclusion
settings. No production process was restarted. Use
`scripts/check_background_control_selection_readback_cases.py` for the corrected
fixture; the initial fixture is preserved because the live checker pins it.

Resource plan: one CPU, 16 GiB RAM, no swap, negligible proof output and an
uncalibrated 0.5–12-hour range. Full production verification, actual match
balance, structural outcome qualification and evolutionary effect estimates
remain pending. No GPU or paid infrastructure is used.

## Full selected-control balance assessment queued

Queued balance assessment behind the verified live full-selection auditor.
All 432 guide/policy/scenario combinations remain, including empty strata.
Eight permutation-invariant pair features are assessed: raw sequence distance,
log positive sequence distance, mean log length, absolute log length ratio,
mean/minimum pLDDT and mean/maximum low-confidence fraction. Positive-log rows
exclude zero-distance pairs without adding an epsilon; raw distances, zero
counts and the per-feature pair denominator remain explicit.

For each feature, report selected target/control means and sample SDs,
pooled-SD standardized mean differences, mean/95th-percentile/maximum paired
absolute differences, and the selected-target mean shift relative to all
modeled target records in that guide, scaled by the original-target SD where
estimable. Missing/one-pair/zero-variance cases receive explicit statuses, not
fabricated standardized values. No automatic balance-pass cutoff is imposed.
Coverage summaries retain matched/unmatched targets, taxa, families, distinct
background nodes, maximum reuse, top-five background concentration and
identical-model target/control counts.

These are event-weighted descriptive diagnostics, not independent-sample tests.
Neither small aggregate differences nor matching itself corrects phylogenetic
or family dependence, structural outcome missingness or prediction uncertainty.
Selection shifts are relative to the observed modeled target universe, not all
526 project taxa. Taxon-balanced and family-aware sensitivity remains necessary.

Known moments/SMDs and selection shifts, empty/single/zero-variance samples,
positive-log exclusions and endpoint swaps passed fixtures. Resources: one CPU,
16 GiB RAM, no swap, 1 GiB output allowance and an uncalibrated 0.1–4 hours after
selection verification. Production balance and independent readback remain
pending. No GPU or paid resources were enabled.

## Full balance verification queued

Queued `fungal-background-control-balance-readback-20260927.service` after
confirming the balance controller's live identity. The independent checker
rejoins all selected target/control identities to frozen node records, rebuilds
all eight feature definitions without importing producer helpers, and uses
separate per-stratum pandas statistics to check every mean, SD, standardized
difference, absolute-difference quantile, baseline shift and nonestimable status.
All 432 coverage cells, 3,456 balance cells, original-target denominators,
selected taxa/families, zero distances, identical models and background reuse
counts are in scope. Empty strata must also be present.

253 independent reconstruction fixtures passed, including empty, single-pair,
constant and missing-value cases. Floating summaries use fixed tight numerical
tolerances; statuses/counts/identities require exact agreement. Resources: one
CPU, 24 GiB RAM, no swap, negligible proof output and an uncalibrated 0.1–4 hours
after the producer. Production balance and verification remain pending; no
biological effect or independent-sample test is claimed. No GPU or paid
infrastructure was enabled.

## Selection and balance completed

Both full independent checks passed. Archived final proofs and artifact hashes
for all 2,786,912 selections, 47,190,168 scenario decisions, 432 coverage summaries
and 3,456 balance summaries. See the [completed balance report](duplication-control-balance-20260927.md)
for illustrative strict-background results, control reuse and the substantial
shift toward more confidently modeled targets. Earlier pending entries above
are stage history. Matched structural effects remain pending.

## Additional coordinate production complete

The additional-background coordinate producer completed all 148,104 models in
149 shards, with every model validated and no content rejections. Its systemd
unit exited successfully. Checked all 149 compressed artifact hashes and their
individual receipts against the final producer receipt; archived the result in
`metadata/background_coordinate_production_completed_20260927.json`.

The independent raw-CIF verifier is now running under its previously queued
process (PID 952285, creation time 1790499798.48). Its final proof remains
pending. Whole-chain input preparation and domain extraction retain their
existing dependency gates; coordinate production alone does not qualify these
structures for evolutionary inference. GPU prediction remains paused.
## Full additional-background coordinate readback passed

The independent coordinate audit completed successfully for all **148,104
additional models and 66,498,349 residues**. Every accepted exported C-alpha
sequence, coordinate and confidence value was reconstructed from the frozen
raw CIF atom rows. The verifier shares the CIF lexical parser with production
but uses a separate extraction/checking implementation. There were no rejected
models in this cohort.

After the audit unit terminated with exit status zero, all 149 individual proof
hashes and compressed coordinate-shard hashes were rechecked, along with pinned
inputs, source bindings and summed model/residue/disposition counts. Completion
is archived in `metadata/background_coordinate_readback_completed_20260927.json`;
the full proof is
`results/structural_comparisons/background-coordinate-readback-20260927-v1/receipt.json`
(SHA256 `9589e8aaedd41a5ef39629301b719041617b561f0cc30c0301dd25055c8b9a7e`).

The existing whole-chain materializer and domain extractor automatically
advanced past their dependency gates. Their exact live process identities were
verified; logs showed 20/149 whole-chain shards prepared and 16,561/131,986
domain intervals processed at the checkpoint. The queued whole-chain alignment
controller also remained live, waiting for its complete input stage. No
duplicate jobs or new structure predictions were launched.

This completes coordinate validation, not structural-distance estimation or
the duplication-effect analysis. Whole-chain and domain comparisons and their
numerical checks remain downstream. GPU prediction remains paused.
## Background comparison inputs complete; alignments running

Whole-chain input production completed all 148,104 additional models under both
masks (296,208 dispositions): 148,104 full inputs and 145,522 pLDDT70 inputs are
ready; 2,582 confidence masks retain too few residues and remain explicit.
The exported PDBs occupy 8,926,165,745 bytes. Verified the complete two-mask
manifest grid, unique identities, sequence/position lengths, all disposition
counts, manifest hash, plan binding and successful producer termination.
This checkpoint does not independently reread every whole-chain PDB; upstream
raw-CIF validation and downstream alignment checks have separate scopes.

Domain input production completed 131,986 intervals from 54,004 models under
both masks (263,972 dispositions). All full intervals and 131,639 pLDDT70
intervals are ready; 347 masked intervals retain too few residues. The domain
PDBs occupy 3,315,586,773 bytes. The independent full domain serialization audit
also passed: every accepted atom identity, original residue position, XYZ,
confidence and occupancy was checked, covering **40,907,133 C-alpha atoms**
across the full and filtered exports. This atom count includes repeated
observations across masks and overlapping intervals, not unique residues.

Completion evidence is archived in
`metadata/background_input_stages_completed_20260927.json` and
`metadata/background_domain_input_completed_readback_20260927.json`.
Producer and domain-auditor units terminated successfully. Both previously
queued comparison controllers advanced automatically and were verified live
by PID, creation time and command. Whole-chain comparisons had processed
448/285,800 directed mask dispositions and domain comparisons 9,344/267,716
at the checkpoint. These counters include explicit unavailable inputs and
are not counts of independently validated evolutionary effects.

Structural-distance execution and numerical readback remain unfinished.
Confidence filtering does not establish PAE reliability or physical domain
boundaries, and full versus masked scores use different length denominators.
No structure predictions or duplicate comparison jobs were started.

### September 28: complete whole-chain background production

The whole-chain producer finished successfully at 21:14 EDT, retaining all
285,800 directed mask/order dispositions for 71,450 new model pairs.
There are 142,900 full-chain alignments, 139,220 pLDDT70 alignments and
3,680 explicitly unavailable pLDDT70 inputs. No native-error, parse-error or
timeout disposition is recorded. These outcomes do not include previously
computed primary/reference pairs, which retain their separate provenance.

The final manifest was checked against the complete expected pair/mask/order
universe, with unique paths, exact disposition counts, manifest checksum and
input-bundle checks. Producer terminal state and receipt are archived in
`metadata/background_alignment_producer_completed_20260928.json`.
This is production completion only. The existing independent auditor has
advanced from waiting to computation and is writing its numerical readback
table. Individual checkpoint provenance, residue mapping, RMSD and sequence
identity checks remain incomplete; no validated whole-chain background
contrast or biological duplication effect is claimed yet.

### Strict audit failure and complete discrepancy census

The strict audit terminated with exit code 1 on September 28 at 21:20 EDT,
after logging 80,000 checked dispositions. It rejected an RMSD difference
larger than the unchanged 0.00501 Å printed-rounding allowance. The first
affected pair is AF-A0A015IZU3-F1 versus AF-A0A397SIQ4-F1 (both version 6),
pLDDT70 mask, order zero. Only two residues align. Native RMSD is 0.00 Å;
coordinate reconstruction gives 0.013327706 Å. The independent two-point
segment-length formula gives the same result. A two-point mapping has a
nonunique rotation; this one case does not characterize all remaining pairs.

The failed unit and partial output are preserved. A separate complete diagnostic
pass now repeats all provenance, mapping and numerical checks, retains every
RMSD discrepancy explicitly, and sets scientific eligibility false. It does
not relax the strict audit or rerun native alignments. Eight-disposition fixtures
passed, including rejection of rehashed corrupt stored metrics; the real first
failure passed the analytic two-point check. The new diagnostic uses one CPU,
16 GiB RAM, no swap, a 1 GiB output allowance and a 0.1–4 hour planning range.
Full discrepancy census and subsequent geometry assessment remain pending.
See `metadata/background_alignment_strict_audit_failure_20260928.json` and
`metadata/background_alignment_rmsd_diagnostic_launch_20260928.json`.

The full geometry assessment is queued behind terminal success of that
diagnostic. It covers all 282,120 successful directed background alignments,
reconstructs their mappings, repeats diagnostic values and measures coordinate
rank and proper-rotation curvature. The established reference-geometry method
is reused in a separate background runner; point, two-point, planar, reflected
and random fixtures passed independent quaternion/alternate-SVD checks, with
corrupt geometry fields rejected. One CPU, 16 GiB RAM, no swap and a 2 GiB
output allowance apply; the broad 0.2–8 hour estimate excludes waiting.
All RMSD flags remain explicit. This assesses numerical identifiability,
not prediction accuracy, alignment homology or biological effects. Full
independent geometry readback is still required. Launch provenance is in
`metadata/background_alignment_geometry_launch_20260928.json`.


### Independent background geometry readback queued (September 28)

The full background geometry output now has a queued independent readback,
`metadata/background_alignment_geometry_readback_plan_20260928.json`.
It waits for the exact geometry producer and requires terminal success before
reading its complete, source-bound output. Every successful alignment is checked
against hashed original coordinates and native mappings. An independent 4x4
quaternion eigensystem checks proper-rotation curvature and optimum; LAPACK
`gesvd` checks coordinate and cross spectra. Complete row membership and counts
are required. Near-zero quaternion gaps are recorded rather than interpreted as
proof of rank at machine precision. RMSD flags remain unchanged.

The new checker passed single-point, two-point, planar, reflected, and general
3D fixtures and rejected altered rank, curvature, classification, and ratio.
The full readback is queued, not completed. Resources are one CPU, 16 GiB RAM,
no swap or GPU, with an estimated 0.2–8 hours after prerequisites finish and up
to 1 GiB output. Launch identity and checksums are recorded in
`metadata/background_alignment_geometry_readback_launch_20260928.json`.


### Full background RMSD census completed and checked (September 28)

The diagnostic service finished successfully at 21:34 EDT. Complete table
membership, plan pins, artifact hashes, counters, all discrepancy fields, and
terminal success were checked by `scripts/check_completed_background_rmsd_diagnostic.py`.
The 285,800 dispositions comprise 282,120 successful alignments and 3,680
explicit unavailable-input dispositions. Of the successful alignments, 282,115
pass the unchanged 0.00501 Å rounding tolerance and five fail. All five are
confidence-masked two-residue mappings across three distinct model pairs.
Maximum RMSD discrepancy is 0.013327706036000187 Å. Across the full table,
19 mappings contain one aligned residue, 58 contain two, and 282,043 contain
at least three.

`scripts/check_background_two_point_discrepancies.py` independently confirmed
all five flagged RMSDs using half the absolute difference between the two CA
segment lengths. This exact two-point formula agrees with the reconstructed
least-squares RMSDs within 1e-10 Å. The native discrepancies remain flagged;
this does not establish their software-level cause. Two-point mappings do not
determine a unique rotation. Full geometry assessment and its independent
readback remain in progress/queued, and scientific eligibility remains false.
Evidence is recorded in `metadata/background_rmsd_diagnostic_completed_20260928.json`
and `metadata/background_two_point_discrepancy_review_20260928.json`.

The completion recorder `scripts/check_completed_background_alignment_geometry.py`
is prepared for the whole-protein geometry audit. It requires terminal success
for both producer and independent audit, verifies their source hashes and plans,
checks every serialized alignment key and RMSD flag against the full diagnostic
table, and checks exact agreement of all degenerate records with the independent
audit. It also records the cross-tabulation of mask, geometry, and RMSD status.
The terminal gate was verified to refuse completion while the producer is active.
No completed geometry record has yet been written; run the recorder only after
both services finish successfully.

### Whole-protein geometry production complete; independent audit running

The producer completed successfully with all 282,120 alignment geometries.
Complete output membership, counters, hashes, and unchanged numeric diagnostic
flags have been checked. It reports 282,043 unique rotations at the numerical
tolerance (142,900 full and 139,143 confidence-masked), and 77 nonunique
confidence-masked rotations. The latter are exactly the 19 one-residue and
58 two-residue mappings. All five RMSD discrepancies belong to this group;
72 other nonunique mappings pass rounding. Rounding agreement alone therefore
does not establish an identifiable rotation.

The independent quaternion audit started automatically and is actively
computing. Numerical uniqueness is not stability under coordinate uncertainty
or biological eligibility. The production-only evidence is recorded in
`metadata/background_alignment_geometry_producer_completed_20260928.json`.
The final completion recorder must wait for the independent audit to finish.

### Isolated local Kabsch source probe

The locally installed aligner directory includes a monolithic `USalign.cpp`.
Its Kabsch routine derives singular values through an analytic eigensystem of
coordinate products and calculates the residual sum of squares from those values.
`scripts/probe_background_kabsch_two_points.py` extracts that exact routine into
a separate C++ harness, leaving the production executable untouched. All five
flagged mappings were tested using original, centered, and translated (+1000 Å)
coordinates, for 15 total cases. The two-point analytic formula was independently
checked for every serialized input/output row. Maximum absolute errors relative
to that formula were 0.0134824 Å (original), 0.0119144 Å (centered), and
0.0133277 Å (translated). Translation-dependent output despite invariant exact
RMSD supports numerical sensitivity in this local routine for rank-deficient fits.

This does not prove the installed executable was built from that source, and
the native pipeline passes already transformed coordinates to its final Kabsch
call. The probe is therefore evidence for a mechanism, not an exact reproduction
or a clearance of the original discrepancies. Compiler version, command,
source/harness/binary checksums and all inputs/outputs are saved under
`results/structural_comparisons/background-kabsch-two-point-probe-20260928-v1`.
The full 15-case check is recorded in
`metadata/background_kabsch_two_point_probe_completed_20260928.json`.

### Background order summary queued

`metadata/background_alignment_usable_orders_plan_20260928.json` now queues the
complete 71,450 newly measured background pairs behind terminal success of the
independent geometry audit. Both masks and orders are retained (142,900 summary
rows). Numerical exclusions require original RMSD rounding agreement, at least
three paired residues, and numerical rotation uniqueness. Excluded metrics stay
blank; native statuses and every applicable reason remain explicit. This is not
coverage/confidence qualification or biological acceptance. Earlier reused
primary/reference pairs still require a separate source-bound join.

The producer normalizes endpoint scores consistently across reversed input order,
then runs an independent complete summary readback checking every metric,
exclusion, blank and order difference. Endpoint normalization and combined
exclusion/blank fixtures passed. Resources: one CPU, 16 GiB RAM, no swap/GPU,
up to 2 GiB output, estimated 0.1–2 hours after prerequisites. The full summary
and its readback have not yet completed. Exact launch identity is saved in
`metadata/background_alignment_usable_orders_launch_20260928.json`.

### Full independent whole-protein geometry audit passed

The independent quaternion/SVD audit completed successfully for all 282,120
background alignments. The final completion recorder passed complete membership,
source binding, geometry category, unchanged RMSD flag, and degenerate-record
checks. Both producer and audit exited successfully. The result confirms
282,043 numerically unique rotations and 77 nonunique rotations, with all five
original RMSD discrepancies preserved. Evidence:
`metadata/background_alignment_geometry_completed_20260928.json` and
`results/structural_comparisons/background-alignment-geometry-readback-20260928-v1.json`.
The queued order-summary stage has started. Coverage, confidence, domain
orientation and biological interpretation still require downstream analyses.

The background order summary and its full independent readback subsequently
completed successfully at 21:59 EDT. All 142,900 pair/mask rows and 3,384,504
numeric values were checked. Full masks retain both usable directions for all
71,450 pairs. Confidence masks retain both directions for 69,570 pairs, one
for three, and neither for 1,877. The latter includes unavailable inputs and
numerical exclusions; it is not a count of RMSD errors. In total, 282,043 directed
alignments pass the numerical screen; 77 aligned directions are excluded (72
short/nonunique and five additionally failing RMSD rounding). Evidence is in
`metadata/background_alignment_usable_orders_completed_20260928.json` and
`metadata/background_alignment_usable_orders_readback_20260928.json`.

### Complete background measurement source union

The audited summaries now cover all 71,461 distinct eligible background pairs:
71,450 new pairs plus the 11 previously measured reference pairs. All 44 reused
checkpoints (two masks × two orders × 11 pairs) were checked against their
original manifest, exact ordered endpoint identities, materialized input files,
executable checksum and options. Original summary fields, including exclusions
and blanks, survive unchanged. The serialized combined table has exactly
142,922 unique pair/mask rows matching the entire inventory; counts and source
assignments were checked. Both full-mask directions pass the numerical screen
for all 11 reused pairs; confidence masks have both directions for four and
neither for seven. Missing directions are preserved without imputation.

Reproduce with `scripts/join_background_alignment_order_sources.py`; output:
`results/structural_comparisons/background-combined-orders-20260928-v1`.
Completion evidence: `metadata/background_combined_orders_completed_20260928.json`.
Candidate/event linkage, identical-model and neither-guide dispositions,
coverage/confidence screens and evolutionary tests remain downstream.

### All background candidates linked to audited measurements

`scripts/link_background_candidate_measurements.py` now preserves all 78,372
background candidates under both masks (156,744 rows), with original gene/taxon
identities, guide membership, native orthology flags, sequence distances,
length/confidence covariates and measurement dispositions unchanged. All 71,461
measured pairs are linked to their audited summaries. The 6,448 identical-model
and 463 neither-guide-eligible candidates retain blank fit fields; structural
identity is not converted to an imputed distance.

Structural metrics retain explicitly named canonical endpoints. Candidate
orientation is recorded as same for 45,704 measured candidates and reversed for
25,757; this prevents accidental assignment of endpoint-normalized coverage or
TM-score to the wrong gene. Every serialized source and measurement field was
checked, followed by independent checks of all 156,744 unique candidate/mask
identities, pair hashes, orientations, mask assignments and unmeasured blanks.

Output: `results/structural_comparisons/background-candidate-measurements-20260928-v1`.
Evidence: `metadata/background_candidate_measurements_completed_20260928.json`.
This preserves sampling denominators; it does not complete coverage/confidence
qualification, phylogenetic matching or evolutionary inference.

### Whole-protein coverage sensitivity screens complete

All 78,372 candidates were screened under both masks using the existing six
whole-protein criteria from `metadata/whole_protein_common_fits_plan_20260927.json`:
30 or 50 aligned residues and 50%, 70% or 90% coverage of both original proteins.
Both input orders must be numerically usable and meet the thresholds. Coverage
of confidence-masked inputs uses full original protein lengths, not retained
lengths. Unmeasured candidates remain explicit failures of measurement eligibility.
No favorable order is selected and no confidence-calibration claim is made.

For at least 50 residues, both masks pass for 45,330 candidates at 50% coverage,
29,430 at 70%, and 9,129 at 90%. These are descriptive sensitivity results,
not independently replicated controls or biological significance. The masks
need not be nested: one candidate passes the 50% masked screen but not the
full screen because the aligner can choose a different mapping after masking.

All 940,464 screen decisions across 156,744 candidate/mask rows were independently
replayed using decimal ceiling cutoffs, including exclusion lists and exact
both-mask intersections. Boundary, reverse-order failure, numerical exclusion,
identical-model, and overlapping-exclusion fixtures passed. Reproduce with
`scripts/screen_background_whole_protein_coverage.py` followed by
`scripts/readback_background_whole_protein_coverage.py`. Output:
`results/structural_comparisons/background-whole-protein-coverage-20260928-v1`;
evidence: `metadata/background_whole_protein_coverage_completed_20260928.json`.
Domain/orientation checks, confidence calibration, phylogenetic matching and
inferential comparisons remain outstanding.


## Expanded catalog refresh (September 30)

The completed September 27 matching/control results retain their older frozen
structure catalog. They are not silently relabeled as expanded results. The
September 30 refresh starts again from **all 140,719 resolved trees and all
3,154,373 guide-specific terminal sister pairs**, using the independently
verified September 28 bridge of 1,955,694 protein/model links. No missing-model
pair is removed from the inventory and no structural response selects a pair.
All prior trees, duplication flags, branch lengths and annotation policies
remain unchanged; only the frozen structure catalog and downstream source
bindings are refreshed.

Eight source-bound stages have launched with captured PID, creation time and
command: full inventory, independent native-tree reconstruction, complete guide
union, independent sorted-cursor union checking, full native orthology scans,
independent reciprocal-record binary searches, model/metadata inventory and its
complete independent reconstruction. The latter stages wait for their recorded
prerequisites and require checksum-bound completion/readback receipts. Their
plans and process identities are listed in
`metadata/expanded_background_refresh_queued_20260930.json`. Successful future
completion must be collected from complete artifacts, proofs and process
journals; an absent transient unit's default success fields alone are insufficient.

The full inventory producer reports 147,159 profile and 147,099 MAFFT cross-taxon
unreported pairs with two distinct models, compared with 71,858/71,822 in the
older inventory. It retains 13,037/13,104 identical-model pairs. These are
**fully independently checked inventory counts**, not ortholog membership,
matched controls, independent biological events or significance. All native
source fields and counts passed the full postorder reconstruction, and exact
producer/checker completion journals and hashes were collected in
`metadata/expanded_terminal_sister_background_inventory_completed_20260930.json`.
The unchanged full inventory denominators are 1,577,204/1,577,169; the two guide
counts overlap and must not be added as independent events.

Full expanded primary duplication annotations already passed independent
registry reconstruction: **276,682 models, 134,812 distinct pairs and 539,248
pair/policy rows**. All four original Pfam annotation alternatives and repeat
occurrences remain explicit. Expanded background annotation production and
its separate raw-registry readback wait for the refreshed model inventory.
For both sets, missing annotations mean unknown; they do not establish domain
absence. Expanded sibling-reference comparisons require a separate refreshed
reference ledger and are not replaced by the primary annotation inventory.
See `metadata/expanded_duplication_domain_controls_completed_20260930.json`.

Sequence/focal-taxon matching-support production and independent vector-mask
readback are queued for all **283,409 modeled duplicate target links** under
three native-background qualification sets (850,227 support rows). Complete
four-policy architecture support and independent reconstruction then cover
3,400,908 target/set/policy rows. Versioned architecture scripts additionally
check every future input against its verified source artifact manifest before
and after calculation, even when the source did not exist when the static plan
was written. A 60-row fixture exercises all five architecture statuses, three
qualification sets and four policies; changed source artifacts and rehashed
false support both fail. Validation evidence is
`metadata/expanded_architecture_source_binding_validation_20260930.json`.

The original 54 fixed matching scenarios, complete unmatched dispositions,
with-replacement dependence and no-rematching-after-screening policy remain
required for expanded matching. The candidate graph and selections must wait
for full architecture support verification; new structural measurements and
all numerical/coverage/confidence/orientation checks remain downstream. Older
model-fitting runs continue on their explicitly frozen inputs. Neither this
refresh nor the older completed matching proves a duplication effect.

Each inventory/union/membership/model stage has a one-CPU/16-GiB/no-swap limit;
annotation and support stages have one CPU/24 GiB/no swap. Per-stage planning
ranges of 0.1–6 hours exclude prerequisite waits and are not calibrated ETAs.
The initial eight-stage output allowance is 30 GiB, with up to 5 GiB for each
annotation or architecture-support producer. New directories preserve all old
and failed runs. The resource plan was recorded before launch in
`metadata/expanded_background_refresh_resources_20260930.json`; later plans
contain their own allowances. No GPU predictions or paid resources launched.

## Expanded support and candidate graph completed September 30

The complete expanded native inventory, guide union, native orthology
membership, model inventories, primary/background annotations and sequence/
architecture support grids have now passed their independent checks.
[Completion evidence](../metadata/expanded_background_matching_support_pipeline_completed_20260930.json)
binds 150 source/artifact hashes and all 16 exact captured-process completion
journals. It covers 140,719 trees, 3,154,373 terminal-pair rows, 1,582,382 union
pairs and **160,415 modeled background candidates**, preserving one-guide,
unreported-parent and native-orthology dispositions. Computational native
orthology and absence of a reported duplication do not establish biological
speciation or independent events.

The background model inventory has **307,693 models**, 305,434 active endpoints
and 303,802 additional models relative to the expanded primary and preserved
old reference inventories. The 146,172 distinct eligible background model
pairs include 146,161 labeled new relative to those two queues and 11 already
in the old reference queue. That label does **not** mean 146,161 native
alignments must be computed: earlier background measurements and the expanded
reference inventory must be checked for overlap and valid reuse first.
All 584,688 background pair/policy annotation records passed separate registry
reconstruction. Missing annotations remain unknown, not domain absence.

Full sequence support covers 850,227 target/qualification rows, and full
architecture support covers 3,400,908 target/qualification/policy rows for
283,409 modeled duplicate target links. From these audited rows the complete
graph contains **318,037 guide-specific background nodes and 5,059,122 edges**,
with all 1,133,636 target/policy dispositions retained. Both possible endpoint
mappings have independently verified length, mean pLDDT and low-confidence
fraction differences, three joint caliper bands, and shared gene/model/sequence
flags. Graph/covariate completion binds 71 source/artifact hashes and all four
captured-process journals in
[the completion record](../metadata/expanded_background_graph_covariates_completed_20260930.json).

The unchanged 54 fixed scenarios require **61,216,344 target/policy/scenario
decisions**. Selection production has finished, with its full independent
ranking, unmatched and reuse-count readback still running at this checkpoint.
Selection uses sequence distance and endpoint metadata; structural outcomes
and later passing screens do not choose or replace controls. Matching balance,
coordinate and numerical quality, confidence/PAE/orientation controls,
phylogenetic dependence and calibrated biological effects remain outstanding.

Full raw-coordinate validation is running over the 303,802 additional models,
133,998,473 residues and 127,840,912,925 coordinate bytes. It uses four CPU
workers, 32 GiB RAM, no swap and 304 thousand-model shards, followed by a
separate independent reader. Its 0.5–24-hour planning range excludes waiting
and is not a calibrated ETA. The resource record predates launch:
`metadata/expanded_background_coordinate_resources_20260930.json`.
Graph/matching stages have one CPU/32 GiB/no swap each. No GPU prediction or
paid resource was launched.

A complete old/new pair-union source screen and independent SQLite checker
have also been added. They compare model/version identity plus both endpoint
catalog coordinate checksums, sequence checksums and lengths, retain every
matching/changed prior source and deduplicate overlaps between primary,
reference and background inventories. The subsequent version additionally
includes the already measured expanded primary queue. These are planning
screens: raw input bytes, residue/mask mapping, executable/settings, order,
checkpoints, numeric exclusions and source proof lineage must pass before any
result is reused. Sources, output locations and launch identities are in
`metadata/full_pair_reuse_plan_20260930_v2.json` and
`metadata/full_pair_reuse_queued_20260930_v2.json`. Neither screen admits or
copies structural results. All eight scientific aims remain incomplete.

### Complete pair-union source inventory verified

The old/new source screen and its extended existing-measurement version both
passed exhaustive independent SQL reconstruction. Completion binds 45
source/artifact hashes and four exact captured-process journals in
[the reuse inventory evidence](../metadata/full_pair_reuse_inventory_completed_20260930.json).
The full current primary/reference/background union contains **336,292
distinct model pairs**. Relative to original collections plus the completed
expanded primary queue, **236,650 have matching catalog source signatures**
and **99,642 are new to those collections**; no source-change-only pairs were
found. These counts are structural comparison pairs, not proteins lacking
predictions, and matching catalog signatures do not authorize result reuse.

| Current role | Distinct pairs | Matching catalog sources | New to existing collections |
|---|---:|---:|---:|
| Expanded primary | 134,812 | 134,812 | 0 |
| Expanded references | 55,701 | 30,754 | 24,947 |
| Expanded backgrounds | 146,172 | 71,460 | 74,712 |
| Distinct union | 336,292 | 236,650 | 99,642 |

Role rows overlap: 17 new reference/background pairs are shared, and existing
pair overlaps are also retained. Never add role counts as independent work or
biological events. Every endpoint source checksum, sequence checksum and length
was compared in numeric model/version order. Actual reuse still requires
raw/materialized coordinate and mask/residue-map equality, executable/options,
input-order/checkpoint bindings, full numerical proofs and preserved exclusions.
The first old-catalog-only inventory is preserved alongside the extended one.

Full additional-reference coordinate validation and separate independent
native-CIF reconstruction are now running for **24,804 models**, 9,982,283
residues and **9,597,401,251 raw bytes**. Each has two CPU workers/8 GiB/no swap;
25 thousand-model shards retain all rejection dispositions. The exact residue
count and maximum length are authoritative in
`metadata/expanded_duplication_reference_coordinate_resources_20260930.json`.
The v2 validator requires the stronger native gene/model ledger readback;
its exact additional-partition and original coordinate algorithms are unchanged.
A full prelaunch partition check passed. Runtime settings for both transient
services were verified against their planned limits; no GPU setting changed.
Coordinates, alignments, common-residue triads and biological asymmetry remain
unqualified. Resource planning predates launch and does not give an ETA.

## Full expanded fixed matching completed September 30

Both selection and independent exhaustive readback have now finished. All
**61,216,344 target/policy/scenario decisions** passed, including **4,250,692
selected records** and **56,965,652 unmatched decisions**. The unchanged design
has 1,133,636 target/policy records and 54 fixed scenarios. Every eligible
candidate, chosen identity/order/score, tie count and runner-up gap was
independently reconstructed, with all per-guide/policy/scenario summaries and
background reuse counts checked. Completion binds 28 source/artifact hashes
and both exact captured-process journals in
[the matching evidence](../metadata/expanded_background_fixed_matching_completed_20260930.json).

These are overlapping scenario decisions and with-replacement selections,
not 4.25 million independent evolutionary events. The design excludes shared
genes/models/sequences, retains exact zero distances without an epsilon, and
uses one endpoint mapping that satisfies all joint cutoffs. No structural
response or later passing screen selected a control. No rematching will occur
after outcome, confidence or coverage filtering.

Full descriptive pre-measurement balance and separate independent reconstruction
have launched for **432 guide/policy/scenario strata and 3,456 feature rows**.
The eight fixed pair-level features cover sequence distance, positive-log
sequence distance, length, length asymmetry, pLDDT and low-confidence fractions.
All original modeled target baselines, matched taxa/families, exact reuse,
zero distances and nonestimable variance outcomes remain explicit. These are
balance/selection diagnostics, not independent-sample tests or biological
effects. The producer has one CPU/16 GiB/no swap; the reader has one CPU/24 GiB/
no swap. Planning ranges are 0.1–8 hours excluding waits, not ETAs. Results are
pending; post-measurement-screen balance is also still required.

Full supplemental alignment input preparation and independent written-input
checks are queued behind complete raw-coordinate readback, retaining all
**303,802 additional background models and 24,804 additional reference models**
and both full/pLDDT70 masks. All short/rejected states stay in their denominators.
The resulting grids have **607,604 background and 49,608 reference input
dispositions**. Overlapping models in the two roles remain separate source
ledgers until exact byte/mask/position equality is verified; these are not
disjoint input collections or independent predictions.

The reference materializer v2 changes only its partition-loader import to
require the stronger native gene/model ledger proof; coordinate/mask serialization
is unchanged. A new full checker independently parses every ready PDB and
verifies original residue position, amino acid, model/version, source checksum,
coordinate/confidence rounding, rejection/short state and all summary counts.
It also checks exact native model partitions, every shard and full proof/hash
lineage before and after reading. Isolated end-to-end fixtures cover both roles,
numeric model versions 6/10 and all six mask/disposition outcomes; rehashed
wrong residue maps, wrong versions, missing inputs and wrong PDB coordinates
fail. These fixtures do not establish production completion or independently
adjudicate raw-CIF rejection causes.

Plans and all four exact queued processes are recorded in
`metadata/expanded_additional_alignment_inputs_queued_20260930.json`.
Background preparation/checking has one CPU/16 GiB each; reference stages have
one CPU/8 GiB each; all have no swap. Resource plans predate launch. Maximal
PDB bytes are bounded by two 81-byte-per-residue masks plus file terminators,
with separate manifest/filesystem allowances: 48 GiB background and 4 GiB
reference producer output allowances. Each stage keeps a 100 GiB disk reserve.
Uncalibrated 0.1–12-hour planning ranges exclude coordinate prerequisite waits.
No native alignment, GPU prediction or paid resource was launched. Actual
native-result reuse/measurement quality, prediction/domain/orientation controls,
phylogenetic dependence and calibrated duplication/asymmetry inference remain
required. All eight scientific aims remain incomplete.

## Expanded pre-measurement balance independently verified

Full balance production and independent direct selection/node joins completed
for **all 4,250,692 selected records, 432 strata and 3,456 eight-feature rows**.
All moments, quantiles, pooled-SD standardized differences, baseline selection
shifts, nonestimable states and coverage/taxon/family/reuse/zero-distance counts
passed. Source/artifact hashes and both exact captured-process completion
journals are collected in
[the completion evidence](../metadata/expanded_background_control_balance_completed_20260930.json).
This is pre-measurement descriptive balance; post-screen balance, predictor and
phylogenetic uncertainty, calibration and effect inference remain required.

Across the full grid, matched counts range from **531 to 26,502 target records**
and matched taxa from **93 to 199**. These are dependent sensitivity strata;
their ranges are not uncertainty intervals or counts of independent events.
The original target pools have 141,724 profile and 141,685 MAFFT records.

![Expanded metadata balance and target selection](figures/expanded_background_control_balance_20260930.png)

The same illustrative S45/alignment-E-value scenario as the older catalog is
shown, without choosing it as a confirmatory primary test. It contains **18,503
matched targets across 198 taxa under each guide**, with 1,935 zero-distance
pairs excluded from positive-log-distance summaries only. Its matched
standardized mean differences range approximately **−0.047 to 0.101**. Relative
to the full modeled target pool, matched targets have mean-pLDDT shifts of
about **0.843 original-target SD** and minimum-pLDDT shifts of about **0.911 SD**;
log-length-asymmetry shifts are about **−0.635 SD**. Therefore close matches
among measured targets do not establish representative sampling of all targets.
The two panels explicitly use different SD denominators; guides are alternatives,
not independent replications. Confidence values are prediction metadata, not
validated accuracy or biological effects.

All 32 figure points and sample counts passed separate source-table readback;
the rendered PNG was visually reviewed, with matching PDF/SVG and full point
[table](figures/expanded_background_control_balance_20260930.tsv) exports.
[Figure evidence](../metadata/expanded_background_control_balance_figure_completed_20260930.json)
records source hashes and visual review. The complete source table remains
`results/orthology/expanded-background-control-balance-20260930-v1/covariate_balance.tsv`;
all strata are retained regardless of illustrative plotting. All eight
scientific aims remain incomplete.


## Full expanded background coordinate producer finished September 30

The raw-coordinate validation producer finished all **303,802 additional
background models / 304 shards**, retaining **303,802 validated dispositions
and no content rejections**. Completion checked **321 source/artifact hashes**,
including every gzip coordinate shard, and the exact captured producer's
successful completion journal:
[producer evidence](../metadata/expanded_background_coordinate_producer_completed_20260930.json).

The independent native-CIF/residue reader is active with its four CPU/32 GiB/no
swap allowance. This producer closure does not imply that the independent
reader or the written-input preparation is complete. The existing full
background input materializer and PDB/residue reader still depend on that
native check. Source/result reuse, overlapping reference/background inputs,
new native background measurements, confidence/coverage/PAE/domain controls
and biological comparison remain required; fixed matches are not rematched
following quality exclusions. No GPU inference or paid resources launched.


## Full supplemental background inputs completed October 1

Every 303,802 additional model disposition and 133,998,473 residue coordinates
passed the full native-CIF/exported-coordinate reader across 304 shards. Every
607,604 full/pLDDT70 written PDB input disposition passed independent checks;
5,752 short pLDDT70 inputs remain explicit. Written bytes total 18,056,589,431.
[Composite completion](../metadata/expanded_background_inputs_completed_20261001.json)
binds 651 source/artifact hashes and all four exact original journals. This
closes input materialization, not biological matching or structural divergence.
The full 146,172-pair background measurement design still needs actual reuse
qualification/new measurements, coverage and predictor/orientation calibration.

The [complete input-partition overlap diagnostic](../metadata/expanded_background_partition_overlap_diagnostic_20261001.json)
exhausted the expanded-primary, old-reference-additional and expanded-background
model partitions. The background partition is disjoint from both others. Six
models occur in both primary and old-reference input collections: all 12
full/pLDDT70 states have identical physical semantic fields, source hashes and
written PDB hashes, checked on all 24 actual files. Collection-specific paths
and coordinate-shards remain separate provenance. The old background handoff
rejects all overlaps, so it cannot be used unchanged for this expanded design.
A new full handoff must permit only explicitly proven byte-equivalent overlaps
and preserve both sources. This narrow exhaustive overlap diagnostic does not
qualify the entire input union or any old alignment/native result reuse.


## Full four-collection input union and new measurements October 1

The frozen design contains 146,172 eligible physical pairs and 305,434 active
models, drawn from 307,693 total catalog models. The active model count is
distinct from the total catalog count. Complete catalog reuse screening finds
74,712 pairs without an old measurement source and 71,460 matching old
catalog pairs. The latter count is a reuse candidate count, not accepted reuse.
The older inventory's `work_disposition` is retained alongside this current
measurement partition; its 146,161 additional-pair flag does not specify how
many genuinely new native alignments are needed.

Three modern input collections alone omit 224 required active models. Preserve
the older reference supplementary collection as a fourth source. Full written
PDB readback passed all 14,540 legacy models/29,080 full and pLDDT70 states:
14,540 full ready, 14,050 masked ready and 490 short masked inputs. It checked
816,376,377 actual written bytes across 15 shards. The
[versioned legacy completion](../metadata/legacy_reference_written_inputs_completed_20261001_v2.json)
binds 71 hashes and both exact original materializer/reader journals. Its first
completion attempt failed because the old materializer receipt had not been
explicitly pinned as evidence. The corrected v2 plan pins both evidence
receipts; the passed input reader was not rerun, and the original failed plan
and launch remain preserved.

The four full input collections contain primary 276,682, modern reference
24,804, additional background 303,802 and legacy reference 14,540 model records,
or 1,239,656 source model/mask records including overlaps. A new immutable
[full input-union plan](../metadata/expanded_background_native_handoff_plan_20261001.json)
checks every source record and globally overlapping semantic projection,
actual PDB bytes, source sequence/coordinate fingerprints and original residue
positions. Selection follows the declared primary/reference/background/legacy
collection order, independent of readiness or result quality, and preserves
all source origins. Modern reference/background overlap includes 1,012 models;
the earlier narrow diagnostic of the older three-collection design does not
contradict this different four-collection scope.

The producer exports every active model, all 610,868 active model/mask states,
all global overlapping input records and all 146,172 work-partition rows. A
complete independent SQLite/source reconstruction and two exact original
completion journals must close the stage. At pipeline launch this input union
was live, not complete. The full union and its source/artifact hash archive
remain outside Git, with small versioned completion locators. The compact
active-input export retains a semantic hash and the exact chosen source row's
path/ordinal/hash rather than copying large residue-position arrays. Subsequent
geometry stages reconstruct all chosen original positions from the immutable
original manifests and verify both full-row and semantic hashes.

The input-union producer subsequently finished, with its exact original journal
completion verified. It reports 14,593 globally overlapping models/29,186 mask
states, 305,434 full ready inputs and 299,632 masked ready plus 5,802 short masked
inputs. The new native partition contains 149,424 full-ready and 144,864 masked-
ready directed states, with 4,560 masked unavailable states. It hashed 657,825
actual ready PDB files, including overlap provenance, and 128,483,810,501 raw
coordinate bytes. These producer results still require the live independent
reader and final closure; they do not qualify old result reuse. The
[runtime checkpoint](../metadata/project_runtime_checkpoint_20261001_v3.json)
records the original completed producer and live reader/waiting pipeline
handles, as well as the preserved likelihood, polynomial, BALiPhy and AFDB jobs.

The [new native pipeline](../metadata/expanded_background_measurement_pipeline_started_20261001.json)
has been launched as exact-dependency wrappers. It starts only after successful
input-union closure and covers all 74,712 genuinely new pairs under both masks
and both input directions: 298,848 dispositions. The full ledger preserves the
other 71,460 pairs as pending actual old-input/checkpoint/numeric/geometry reuse
checks. It does not promote a catalog match to an accepted measurement.

The native runner uses the pinned existing USalign binary and settings
`-mol prot -mm 0 -outfmt 0 -ter 2`, saves raw native output and timing, verifies
actual input hashes around each call, and checks saved checkpoint bindings on
replay. There is no automatic retry or source/coverage substitution. Every
short input, native failure, parse error and timeout stays in the full grid.
Numerical/geometry assessment then exhausts every checkpoint, reconstructs
aligned sequence mappings, original residues, identity, confidence, coverage,
proper-rotation least-squares RMSD, coordinate rank and rotation curvature.
RMSD differences beyond the native printed rounding, fewer than three aligned
pairs and nonunique rotations remain explicit numerical exclusions.

The separate reader rebuilds the mappings using vectorized residue indexes,
checks native numeric text separately, uses quaternion rotations for RMSD and
LAPACK `gesvd`/a quaternion eigensystem for spectra and curvature, and checks
every nullable field, disposition, model role, exclusion and summary count.
Near-zero quaternion gaps are recorded without claiming precision sufficient
to resolve ambiguous ranks. TM-scores bind the native text; their optimization
is not independently repeated. Fused numerical/geometry assessment is one
service, so final measurement closure requires three exact original journals:
native, assessment and independent reader. Neither a completed native receipt
nor a default inactive/success/0 systemd observation alone qualifies the work.

The [software validation](../metadata/expanded_background_measurement_fixture_validation_20261001.json)
used five synthetic 20-residue models across four collections, including
overlaps, a legacy-only model, a short masked input and a collinear fit. Actual
USalign computed 12 directed states, of which 10 aligned; the pending synthetic
catalog pair was never computed. Quaternion RMSDs agreed within 2.178e-15 Å;
scaled curvature difference was at most 1.705e-16. Twelve rehashed false exports
were rejected. Checkpoint replay preserved bytes; synthetic native failure,
parse failure, timeout, RMSD discrepancy and degenerate geometry states stayed
explicit. Original-position source corruption was rejected. Only closure
provenance was synthesized in this software suite; it does not test real
completion journals, establish full-data source qualification or constitute
a biological pilot.

[Resources were estimated before launch](../metadata/expanded_background_measurement_resources_20261001.json):
eight native CPU workers/32 GiB, serial two-CPU/64 GiB assessment, reader and
completion stages, no swap, one BLAS thread, 96 GiB output allowance and 100 GiB
disk reserve. The older 285,800-disposition whole-background run used
74.7256 CPU hours. Scaling only by disposition count gives 78.1372 CPU hours
or 9.7671 ideal wall hours on eight workers. This is not an ETA: new lengths,
input availability, I/O and shared CPU use differ, and complete source hashing
adds overhead. The provisional native planning range is 10–72 hours, with
0.5–24 hours per audit stage; it is uncalibrated. A 600-second timeout applies
to each ready native call, not the complete run. The extreme all-calls-time-out
budget is 6,226 wall hours at eight workers, also not a prediction.

Production input union, new measurements and real old-result reuse qualification
are pending. Full confidence/original-length coverage, prediction/PAE/domain
orientation and sequence-locked controls, missingness/phylogenetic adjustment,
matched comparisons and calibrated evolutionary effects remain required.
All eight scientific aims remain incomplete. Existing likelihood, polynomial,
BALiPhy and AFDB jobs are preserved; GPU inference remains paused and there
are no new charges.


## Full input union completed and actual old-result reuse October 1

The entire independent source/SQLite reconstruction passed every 1,239,656
source input state, all 610,868 active states, 14,593 globally overlapping
models and all 146,172 background pairs. The
[completed input union](../metadata/expanded_background_native_handoff_completed_20261001.json)
binds 964,012 actual source/artifact hashes and both original producer/reader
journals. This verifies all 305,434 active raw coordinate models and 657,825
written PDB files including overlap provenance, covering 128,483,810,501 raw
coordinate bytes and 19,651,660,401 PDB bytes. Full/masked ready directed native
states number 149,424/144,864, with 4,560 unavailable masked dispositions. The
eight-worker new alignment runner has started its complete source-input checks;
native measurement, numerical/geometry readback and measurement closure remain
pending. Finishing this input gate does not establish a comparative effect.

The [full reuse pipeline](../metadata/expanded_background_reuse_pipeline_started_20261001.json)
is also launched. It exhausts 146,172 pairs under both masks and directions,
or 584,688 states. A fixed reference-before-background source preference assigns
11 old-reference and 71,449 old-background pairs, or 285,840 candidate states;
the other 298,848 states retain new-native-pending status in this separate
ledger. Both source choices and their original order are explicit. Alternative
matching source labels remain in the full output; a reference input mismatch
does not trigger a fallback chosen after observing an outcome.

The old-background loader validates all three original input partitions and
the complete original native manifest, binds the precise original input bundle,
and checks the full diagnostic/geometry/independent-reader lineage. The
old-reference source uses the existing audited source contract. Every needed
current/old raw source, original residue map, full and retained sequence,
length, mask/status/reason, written PDB bytes, native binary/options/timeout and
directed checkpoint is compared. The candidate pairs use 142,906 current
models (22 reference and 142,884 background), representing 64,779,073 residues
and 61,664,784,336 raw bytes by prelaunch file statistics. Statistics are a
resource estimate; production qualification requires actual byte hashes.

For every compatible reused alignment, the producer reconstructs the original
numeric values on the current byte-identical coordinates using proper-rotation
SVD fits. The reader independently rebuilds every input-identity and directed
reuse export, fixed source choice and old-to-current endpoint order. Vectorized
residue mapping, quaternion RMSD, LAPACK `gesvd` spectra and quaternion
curvature check each reused successful original result. Both full and masked
states remain explicit; confidence uses rounded PDB values. Original RMSD
discrepancies, short alignments, nonunique rotations, unavailable inputs and
native/parse/timeout failures remain in the results. No source substitution,
fresh native optimization or automatic retry is performed. TM-scores bind
native text without independent optimization.

Incompatible input states receive `incompatible_inputs_require_new_measurement`
and remain unqualified, with null result fields rather than zero distances.
If the production check finds any such states, they require a separate explicit
native workload and cannot be resolved by editing the already launched plans.
All original catalog and work dispositions stay in the copied complete ledger.
The reuse producer/reader completion is not the union of reused and new results;
that full union and its independent readback remain required.

The [original source-journal audit](../metadata/expanded_background_reuse_original_source_journals_20261001.json)
verified all eight original reference/background native, numeric, geometry and
independent-reader process invocations, their pinned source plans, original
PID/command journal messages, matching completion resources and successful
terminal states. This journal-only check does not establish current coordinate
or checkpoint identity. Final reuse closure requires those eight journals plus
both new producer/reader journals and the complete current byte/artifact proofs:
ten original journals in total.

The [new software suite](../metadata/expanded_background_reuse_fixture_validation_20261001.json)
passed a full 16-state synthetic grid: 12 old-source states and four new-pending
states, with eight original native alignments reconstructed. It exercises the
actual old-background three-partition I/O contract and reference order reversal.
Six nonunique fits, one RMSD discrepancy, four unavailable states and one usable
state remain explicit. Twelve rehashed false exports were rejected. A separate
mask incompatibility left four states requiring fresh measurement, eight reused
and four new-pending without alternative-source selection. Independent RMSDs
agreed within 2.178e-15 Å; scaled curvature error was at most 1.705e-16. Current
input-union I/O and original proof/journal provenance are synthetic only in this
software suite. It does not establish real production reuse, actual original
journal completion or a biological pilot.

[Prelaunch resource estimates](../metadata/expanded_background_reuse_resources_20261001.json)
allocate two CPU/64 GiB per serial producer, reader and completion stage, no
swap, one BLAS thread, 32 GiB output allowance and 100 GiB disk reserve. The
per-stage planning range is 1–24 hours and is uncalibrated, not an ETA. Full
current-source hashing, original-position reconstruction, selected old-source
readback and current-coordinate fits add I/O and CPU overhead. Existing
likelihood/polynomial/BALiPhy/AFDB handles are preserved; native new measurements
and old-result qualification have separate immutable outputs. GPU inference
remains paused; no new charges. Production reuse and measurement closure,
full union, coverage/domain/PAE/prediction/sequence-locked/missingness/phylogenetic
controls and calibrated effects remain required. All eight aims remain incomplete.

The [runtime check](../metadata/project_runtime_checkpoint_20261001_v4.json)
binds the closed input/archive/proof metadata and records the exact seven live
native/reuse pipeline wrapper handles and five preserved original jobs. Native
USalign execution subsequently began. A
[live progress snapshot](../metadata/expanded_background_native_progress_snapshot_20261001.json)
records saved native dispositions and retains the live wrapper identities.
These changing checkpoint counts are descriptive progress, not full numerical
acceptance. The complete geometry/readback/closure requirements remain in force.


## Full background measurement union queued — October 1

The old-result reuse producer finished all 146,172 physical pairs and 584,688
mask/order states. It classified 285,840 states as matching actual old inputs
and retained results, with 298,848 new-native states pending. No incompatible
input states were reported. The producer reports 282,069 usable old states,
3,694 unavailable states and 77 numerical exclusions, including five RMSD
discrepancies. These are producer counts pending complete independent
quaternion/alternate-SVD readback and ten-original-journal closure. New native
comparisons are advancing under the original eight-worker run. The latest
[exact process and progress observation](../metadata/project_runtime_checkpoint_20261001_v5.json)
preserves the five original likelihood/BALiPhy/AFDB jobs and distinguishes
transient counters from completed source/artifact validation.

The [full union plan](../metadata/background_measurement_union_plan_20261001.json)
covers every 146,172 pair and 584,688 state. It waits for both independently
closed sources: 74,712 new pairs/298,848 states with three original native,
numeric/geometry and reader journals; 71,460 reused pairs/285,840 states with
ten original journals. Fixed old-reference/old-background ownership remains
11/71,449 pairs. Incompatible reuse states block the union until additional
measurements have an explicit closed workload; there is no alternative-source
selection, automatic retry, favorable filtering or missing-as-zero substitution.

The producer streams every full reuse-ledger state and chooses its previously
fixed old source or corresponding new measurement. Its compressed export
retains actual checkpoint path/hash, original source order, directed model
roles, original numeric/geometry rows and flags, native metrics/elapsed time/
return code/error/timeout, plus explicitly typed numeric/geometry fields. An
independent reader indexes every new measurement in SQLite and reconstructs
every export field, complete grid, count and source binding from originals.
The union stage does not repeat native optimization or geometry estimation;
those belong to the completed source gates. Two additional exact original
producer/reader journals and full artifact/source bytes are required for union
closure. All large exports, complete hash dictionaries and archives remain
outside Git under
`results/structural_comparisons/full-background-measurement-union-20261001-v1`.
The future small completion locator is
`metadata/full_background_measurement_union_completed_20261001.json`; its
absence currently means completion is pending.

The [software fixture receipt](../metadata/background_measurement_union_fixture_validation_20261001.json)
passed all 12 states of three synthetic pairs, including reversed old-reference
source orders, both masks, original string metrics versus normalized values,
RMSD/nonunique/short exclusions and unavailable/native/parse/timeout states.
All 15 coherently rehashed false exports were rejected: missing/duplicate rows,
wrong source/direction/model roles/checkpoint, changed original or normalized
metrics/geometry, cleared exclusions, promoted unavailable values, altered
native failures and changed totals. A coherently rehashed incompatible-source
closure failed before output creation. The original records and prior journal
contracts in this fixture are explicitly synthetic: it tests software/source
I/O and export fidelity, not physical measurement correctness, production
completion, original journal verification or a biological pilot. Rerun with
`python scripts/check_background_measurement_union_cases.py`.

[Resources were recorded before launch](../metadata/background_measurement_union_resources_20261001.json):
two CPU/64 GiB per serial stage, no swap, one BLAS thread, 32 GiB output/scratch
allowance and a 100 GiB disk reserve. Full prior hash archives may exceed a
million bindings; every source byte is rechecked. The 0.5–24 hour interval per
stage is uncalibrated planning, not an ETA. The
[queued pipeline](../metadata/background_measurement_union_pipeline_queued_20261001.json)
records the frozen commands/plans, exact original PID/create/CMD captures and
actual cgroup limits. Existing runs were not restarted. No GPU or paid resource
is used; protein structure prediction remains paused.

Production new/reuse/union completion remains pending. The merged measurement
ledger will still require fixed-match coverage/confidence/missingness screening,
sequence-locked and common-residue comparisons, domain/orientation/PAE/predictor
controls, gene/species phylogenetic qualification and comparative calibration.
Numerical availability is not evidence of a duplication effect, branch
acceleration, selection or functional change. All eight scientific aims remain
incomplete.


## Complete background coverage screen queued — October 1

Full independent old-result reuse readback passed all 584,688 design states,
including actual-current-input/checkpoint identity, original numerical flags,
282,146 aligned numerical reconstructions and independent quaternion/alternate-
SVD geometry. The 285,840 retained old states include 282,069 usable, 3,694
unavailable and 77 numerical exclusions, with five RMSD discrepancies. New
298,848 states remain pending native measurement. The ten-original-journal/full-
hash reuse closure is running; successful readback is not itself that closure.
The [current exact runtime evidence](../metadata/project_runtime_checkpoint_20261001_v6.json)
binds the reader receipt and its exact original journal while preserving
original likelihood/BALiPhy/AFDB processes and ongoing new-native measurements.

The [complete background coverage plan](../metadata/full_background_coverage_plan_20261001.json)
uses every 146,172 background physical pair, 292,344 pair/mask row and 584,688
original measurement state. Its six unchanged target screens require 30 or
50 aligned residues and 50%, 70% or 90% coverage of both original full proteins.
Both native directions must be numerically usable and pass; pLDDT70 retained
lengths never replace original protein lengths. Short/coverage reasons overlap,
and unusable orders retain their original native/numerical exclusions.
The producer exports every pair/mask row, original source/order/checkpoint
identity, endpoint versions, full lengths/confidence metadata and both native
order states, including nullable lengths/coverage. Both-mask intersections
remain separate sensitivity counts.

Prelaunch full-ledger enumeration distinguishes 292,326 background endpoint
models from the broader 305,434-model input inventory shared across comparison
types. The complete inventory remains bound; its additional models are not
misclassified as missing background measurements. This distinction was fixed
and tested before first production launch. The current software fixture includes
an unused inventory model, as well as an exact 70/101 versus 71/101 boundary,
confidence masking, asymmetric directions, RMSD/nonunique flags and unavailable/
native-error/timeout states.

The independent reader indexes all original union states in SQLite and rebuilds
every field, direction, exclusion, six-screen decision and both-mask intersection
using decimal ceiling cutoffs rather than producer fraction comparisons. All
1,754,064 screen decisions must be checked. A 36-row physical count table joins
the already closed 134,812-target-pair screen summaries with background full,
pLDDT70 and both-mask counts; all its identities, numerators and denominators
are independently checked. These classes have different sampling denominators;
the table is descriptive, not a matched duplication or ecological effect.

[Software checks](../metadata/full_background_coverage_fixture_validation_20261001_v2.json)
passed all 12 directed states/six pair-mask rows/36 decisions. Thirteen rehashed
false exports were rejected, including altered original or masked denominators,
a boundary pass, cleared numerical/native flags, changed roles/source order/
checkpoint, missing/duplicate rows, changed counts and a changed target-count
comparison. Prior measurement, target and journal contracts are explicitly
synthetic fixture data, not production physics, journal verification or a pilot.
The earlier prelaunch fixture is archived outside Git under
`results/software-checks/full-background-coverage-prelaunch-20261001-v1.json`;
the versioned current receipt refers to the corrected code. Reproduce with
`python scripts/check_full_background_coverage_cases.py`.

The [three-stage pipeline](../metadata/full_background_coverage_pipeline_queued_20261001.json)
waits for the original full measurement-union closure, then runs coverage,
independent readback and complete source/artifact/two-original-journal closure.
All plans/code are frozen after launch; exact original PID/create/CMD and actual
cgroup limits are recorded. [Prelaunch resources](../metadata/full_background_coverage_resources_20261001.json)
allocate two CPU/64 GiB per serial stage, no swap, one BLAS thread, a 32 GiB
output/SQL/hash-archive allowance and 100 GiB disk reserve. Full prior source
hashes may exceed a million bindings and are rechecked. The 0.5–24 hour
planning interval per stage is uncalibrated, not an ETA. Native inference,
GPU resources and new charges are not used. Production outputs stay outside
Git under `results/structural_comparisons/full-background-coverage-20261001-v1`;
the future small completion locator is
`metadata/full_background_coverage_completed_20261001.json`.

Production coverage is pending. Integration must retain all 4,250,692 frozen
selected and 56,965,652 unmatched scenario decisions across 1,133,636 target/
policy records and 54 scenarios; no rematching after numerical or coverage
filtering. Selection reuse, target/control joint attrition, denominators and
post-screen covariate balance remain required, as do common-residue/sequence-
locked/domain/orientation/PAE/predictor/missingness controls, phylogenetic
qualification and calibrated effects. All eight scientific aims remain
incomplete; existing CPU jobs continue and GPU protein prediction stays paused.


## Complete fixed matching coverage queued — October 1

The [old-result reuse closure](../metadata/expanded_background_reuse_completed_20261001.json)
is complete: all 584,688 full-design states, 282,146 numerical alignment
reconstructions, 1,531,355 source bindings and ten original completion journals.
The closure retains 285,840 qualified old states and 298,848 new-native states
pending measurement. Native/numerical exclusions remain dispositions, not
zero structural responses. The full hash archive remains outside Git.

A subsequent [journal audit](../metadata/expanded_background_reuse_journal_revalidation_20261001.json)
rechecked exact original PID/create/command identity, terminal absence and all
recorded completion resource fields for every original service. Raw journal JSON
serialization hashes differed for all ten queries; original hashes and archives
were preserved. Two fresh queries of one completed unit had identical parsed
messages but different bytes ([diagnostic](../metadata/journal_serialization_diagnostic_20261001.json)).
The audit validates original completion fields; it does not claim full canonical
equivalence to the unavailable original raw journal bytes or repeat the prior
complete scientific-input hash/geometry audit.

The [complete node census](../metadata/full_matched_coverage_input_census_verified_20261001_v3.json)
passed every 283,409 target and 318,037 background node, all original catalogue
fields/model versions/gene-to-model roles, all 134,812 target physical pairs and
146,172 background physical pairs. Both complete original model catalogues
contain 584,276 unique identities; their 276,682 target and 307,693 background
rows overlap at 99 identities. All 305,434 current background-input inventory
models were read separately. Every background endpoint is present there; only
146 target endpoint occurrences overlap that inventory. Target-only catalogue
entries are checked against their original catalogue and target physical table.
Every same-model node is retained and catalogue-checked without inventing an
alignment. No confidence-descriptor differences were found. Complete census
file hashes and its exact original v3 process completion are bound.

Two failed census versions remain immutable. Version 1 requested the graph
receipt from the fixed-matching closure, which records selected matching
exports. The dedicated graph/covariate closure provides that receipt and four
original journals. Version 2 correctly used this source but incorrectly required
target-only models in the background inventory. Version 3 checks each complete
source universe in its actual role and all physical membership. These failures
did not stop or restart the original alignments, likelihood or retrieval jobs.

The [integration plan](../metadata/full_matched_coverage_plan_20261001.json)
preserves every original chosen-control row, original endpoint correspondence,
score, ties and all matched/unmatched scenario lists. It covers 1,133,636
target/policy records, 54 scenarios and all 61,216,344 scenario decisions:
4,250,692 selected and 56,965,652 unmatched. Coverage flags are joined by model
ID/version to both full and pLDDT70 physical comparisons; the both-mask screen
is their intersection. Same-model nodes retain an explicit no-alignment
disposition and fail physical eligibility rather than becoming zero effects.
No targets or chosen controls are removed or rematched after seeing outcomes.

Full exports comprise `selected_pair_coverage.tsv.gz`,
`target_policy_coverage_status.tsv.gz` and `matched_attrition_counts.tsv`,
under `results/structural_comparisons/full-fixed-matched-coverage-20261001-v1`.
The complete 7,776 guide/policy/scenario/mask/screen cells include empty strata,
all targets, matched/unmatched counts, target pass counts in all three groups,
control and joint passes, target-only/control-only/neither-pass categories.
Compact six-screen bit fields represent all 76,512,456 selected screen cells;
retained scenario lists represent 1,101,894,192 logical full-design screen cells
without materializing billions of redundant rows. These are dependent design
records, not independent evolutionary replicates.

The reader preserves and compares every original export field, independently
projects model roles/masks in SQLite, verifies every frozen membership and
reconstructs all summary cells. [Software cases](../metadata/full_matched_coverage_fixture_validation_20261001_v3.json)
passed all 54 scenarios, role reversal/version 6/10, mask intersections,
same-model exclusions, all four attrition categories and empty strata; 17
rehashed false exports were rejected. [Separate proof regression cases](../metadata/full_matched_graph_fixture_validation_20261001.json)
rejected eight wrong graph/covariate/matching contracts, including using the
matching closure as graph proof. Prior physics, census and journals are
explicitly synthetic contracts in these software fixtures. Production proof
requires actual closed sources and original journals; these tests are not
physical validation or a biological pilot.

Earlier prelaunch software receipts are preserved outside Git at
`results/software-checks/full-matched-coverage-prelaunch-20261001-v1.json`
and `full-matched-coverage-prelaunch-20261001-v2.json`. The versioned v3
receipt and pinned current scripts include the closed whole-node census gate.

The [queued pipeline](../metadata/full_matched_coverage_pipeline_queued_20261001.json)
waits for complete background coverage closure before integration, then runs
independent readback and source/artifact/two-original-journal closure. Its
[resources](../metadata/full_matched_coverage_resources_20261001.json) were
estimated before launch: two CPU equivalents/64 GiB per serial stage, no swap,
one BLAS thread, 64 GiB export/SQL scratch allowance and 100 GiB disk reserve.
The 1–24 hour planning interval per stage is uncalibrated, not an ETA. Exact
original PID/create/command and actual cgroup limits are recorded; launched
scripts and plans are frozen. No GPU prediction or paid resources were started.
The future completion locator is
`metadata/full_matched_coverage_completed_20261001.json`.

This integration remains pending and does not complete duplication inference.
Post-screen feature balance, background reuse/weights, shared ancestry and
family/taxon/gene dependence, predictor/circularity/PAE/domain/orientation
controls, missingness sensitivity and calibrated effects remain required.
All eight scientific aims remain incomplete.
