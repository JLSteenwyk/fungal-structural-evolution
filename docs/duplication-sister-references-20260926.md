# Provisional sister-clade references for duplicate-pair comparisons

This full-cohort inventory supports the planned structural-asymmetry work.
Both guides' complete reviewed terminal-candidate sets are included: 109,245
profile and 109,228 MAFFT event rows. It does not alter the running frozen
duplicate-pair alignment queue or infer new structures.

For each exact duplicate tip pair, the procedure inspects the immediate
parent and all sister-clade descendants. Every candidate retains sister-gene
and taxon counts, focal-taxon presence, model availability, parent branching
and whether that parent appears in the native duplication table. A
provisional reference is available only if the parent is bifurcating, is not
reported as a duplication, the sister clade has no focal-taxon genes, and
at least one nonfocal sister protein has a frozen structure model.

Missing parents, reported parent duplications, multifurcations, focal-taxon
sister genes and absent models remain explicit dispositions. The absence of
a reported parent duplication does not establish speciation. These are
extant reference candidates, not validated orthologous outgroups or
reconstructed ancestral states. Hidden paralogy, differential loss, uncertain
rooting and topology still need review.

The nearest modeled nonfocal sister tip is identified using the saved
sequence-tree path length from the duplicate node. All ties within absolute
1e-12 are retained; the lexical first tied gene is a deterministic provisional
representative. Selection uses sequence-tree distance and model availability,
not observed structural divergence. Path lengths to both duplicate tips and
the duplicate-pair path are retained as relative sequence divergence, not
calendar time. References are recorded even when another eligibility flag
prevents provisional use, so these cases remain inspectable.

Implementation: `scripts/inventory_duplication_sister_references.py`.
Tests: `python scripts/check_duplication_sister_references.py`. Known-tree
checks passed expected reference eligibility and four exact path distances;
tied tips, zero branch lengths, missing parent, multifurcation, reported
parent duplication and focal-taxon sister cases were checked. Negative
branches and duplicate labels were rejected.

Plan: `metadata/duplication_sister_reference_plan_20260926.json`.
Service: `fungal-duplication-sister-references-20260926.service`.
Launch identity: `metadata/duplication_sister_reference_launch_20260926.json`.
Output: `results/orthology/duplication-sister-reference-inventory-20260926-v1/`.
The full scan uses one CPU, 8 GiB RAM and no swap, with a 1-GiB output
allowance and uncalibrated 0.1–8-hour planning range. Source event tables,
resolved trees, model bridge and completed candidate review are pinned and
rechecked. All candidate rows must receive dispositions before completion.

The inventory is running. Independent readback, cross-guide reference
agreement, reference model quality and domain checks, structural comparisons
to references, alternative-reference sensitivity and asymmetry tests remain
pending. This step does not establish any duplication-associated effect.


## Completed inventory and guide comparison

Both full scans completed and the service exited successfully. Every candidate
received a disposition:

| Immediate-reference disposition | Profile | MAFFT |
| --- | ---: | ---: |
| Provisional reference available | 17,461 | 17,448 |
| Parent reported as a duplication | 54,891 | 54,884 |
| No modeled nonfocal sister protein | 36,817 | 36,820 |
| Parent not bifurcating | 76 | 76 |
| Total candidate rows | 109,245 | 109,228 |

The categories follow the documented priority order; for example, an event
with a reported parent duplication may also lack a modeled reference. They
are not mutually exclusive biological explanations. The table reports the
single output disposition for each candidate.

`scripts/compare_duplication_sister_references.py` independently checked all
original candidate columns and the complete protein-pair universe against
the prior reviewed exports. Among 108,918 shared duplicate pairs, 17,392 are
provisionally eligible under both guides. For 17,365, the complete nearest
reference sets agree and the chosen reference gene/model also agrees. The
remaining 27 have disjoint nearest-reference sets; none have partially
overlapping sets. These cases remain explicit for later reference sensitivity.

Among shared pairs, 28 are provisionally eligible only in profile (13 lack a
modeled reference in MAFFT; 15 have a reported parent duplication there).
Nineteen are provisionally eligible only in MAFFT (one lacks a modeled
reference in profile; 18 have a reported parent duplication there). Separately,
41 profile-only and 37 MAFFT-only duplicate candidates have provisional
references. The complete eligibility matrix is
[`metadata/duplication_reference_guide_eligibility_matrix_20260926.tsv`](../metadata/duplication_reference_guide_eligibility_matrix_20260926.tsv).

Source-plan pins, receipts, all output artifact hashes and count identities
were rechecked and archived in
`metadata/duplication_sister_reference_completed_20260926.json`. This comparison
does not independently reconstruct sister-clade selection or path distances.
High cross-guide agreement is not proof of orthologous reference assignment
or a structural-asymmetry effect. Most candidate pairs still need deeper
gene-tree/context review or additional reference coverage; they remain in
the broader duplication project.

Reference-set fixtures passed shared/guide-only, unavailable, identical,
overlapping and disjoint cases; an empty provisional set was rejected.
Reproduce with `python scripts/check_duplication_reference_comparison.py`.
Reproduce the full comparison with:

```bash
python scripts/compare_duplication_sister_references.py \
  --inventory results/orthology/duplication-sister-reference-inventory-20260926-v1 \
  --plan metadata/duplication_sister_reference_plan_20260926.json \
  --review results/orthology/duplication-candidate-tree-review-20260926-v1 \
  --output <fresh-output-directory>
```

Current comparison output is
`results/orthology/duplication-reference-guide-comparison-20260926-v1/`.


## Complete tied-reference structural workload

`prepare_duplication_reference_comparisons.py` expands every provisionally
eligible event to both duplicate copies against every tied nearest reference.
The 17,461 profile events contribute 18,478 event/reference links; 17,448
MAFFT events contribute 18,466. Expanding both duplicate copies yields
73,888 guide-specific event/reference/copy records. Lexical representatives
are flagged but do not exclude other tied references.

The ledger contains 649 identical-model records, 641 records referring to
model pairs already in the running duplicate-pair queue, and 72,598 records
requiring additional pair computations. After model/version deduplication,
there are 32,862 distinct reference-comparison pairs: 321 already covered by
the frozen queue and 32,541 additional pairs. Full and pLDDT>=70 comparisons
in both input orders would add 130,164 directed dispositions. These are
computational counts, not independent biological observations.

The reference-comparison set uses 49,334 unique models. Of these, 14,540 are
outside the current coordinate queue, requiring validation of an additional
5,784,554,689 bytes (5.4 GiB) of raw coordinates. The existing running queue
and its pins remain unchanged. These extra models have only catalog-level
provenance checks so far; raw coordinates, reference orthology and biological
asymmetry are not validated by this inventory.

`readback_duplication_reference_comparisons.py` independently reconstructed
the complete event × tied-reference × duplicate-copy ledger from the source
reference tables and verified all 73,888 records. It also checked pair keys,
the exact model/pair sets, existing-versus-additional work partitions and
artifact hashes. The producer checked each focal/reference model/version,
sequence hash and path against the frozen bridge/catalog. The independent
readback does not redo that source gene-to-model join or the sister-reference
selection itself. Both receipts are archived in
`metadata/duplication_reference_comparison_completed_20260926.json`.

Output is
`results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1/`,
including `event_reference_comparisons.tsv`, `model_pairs.tsv`, `models.jsonl`
and `additional_models.jsonl`. Reproduce the inventory with:

```bash
python scripts/prepare_duplication_reference_comparisons.py \
  --inventory results/orthology/duplication-sister-reference-inventory-20260926-v1 \
  --inventory-plan metadata/duplication_sister_reference_plan_20260926.json \
  --base-queue results/structural_comparisons/duplication-model-pair-queue-20260926-v1 \
  --catalog results/structures/whole-proteome-afdb-catalog-20260922-v1 \
  --output <fresh-output-directory>
```

Run the readback with `scripts/readback_duplication_reference_comparisons.py`,
passing `--inventory` as the new comparison inventory, `--references` as the
sister-reference inventory, `--base-queue` as the unchanged duplicate-pair
queue, and `--output` as a fresh JSON proof path.


## Additional coordinate validation launched

The supplementary driver `scripts/validate_duplication_reference_coordinates.py`
now validates all 14,540 additional models in 15 checkpointed shards. It
requires the successful inventory ledger readback and verifies the exact
model-record set difference against the original queue, including unchanged
provenance for reused models. It calls the same raw-coordinate validation
function as the primary run without modifying primary scripts or plans.

The resource plan allocates two CPU workers, 8 GiB RAM, no swap, a 4 GiB
estimated output allowance and a 100 GiB free-disk reserve. Maximum sequence
length is 2,271 residues. The uncalibrated planning interval is 0.1–8 hours;
this is not a measured completion estimate. GPUs and paid services are not
used. Validation is running, not yet complete.

Reproduce with:

```bash
python scripts/validate_duplication_reference_coordinates.py \
  --plan metadata/duplication_reference_coordinate_validation_plan_20260926.json
```

The corresponding launch record saves the service, PID, process creation time,
command and plan hash. Output is
`results/structural_comparisons/duplication-reference-coordinate-validation-20260926-v1/`.

The existing independent raw-CIF readback script is queued behind this exact
producer identity using
`metadata/duplication_reference_coordinate_readback_plan_20260926.json`.
It requires the matching successful producer receipt and checks every accepted
exported C-alpha sequence, coordinate and confidence value from source CIFs.
The readback has its own two-CPU/8-GiB limit and writes to
`results/structural_comparisons/duplication-reference-coordinate-readback-20260926-v1/`.
Rejected-model identities/reasons are checked but rejection causes are not
independently adjudicated. Both stages must finish before alignment inputs
can be approved.

The supplementary fixture passed an actual-coordinate producer/readback
handoff, checkpoint reuse, rejection of an incorrectly partitioned model
despite updated artifact hashes, and rejection of an audit bound to another
receipt. Reproduce with:

```bash
python scripts/check_duplication_reference_coordinates.py \
  --models results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1/additional_models.jsonl
```

Reference alignments, independent sister-choice/path validation and biological
asymmetry tests remain pending.

## Independent reference-choice and path readback completed

`scripts/readback_duplication_sister_references.py` reconstructs every output
field for all 218,473 candidates using leaf-interval subtraction for sister
membership and upward edge sums for reference distances. This is separate
from the producer's downward sister traversal. It shares the Bio.Phylo
Newick parser, source trees and native duplication calls. It independently
joins reference model/version identities from the frozen bridge using taxon
and protein names, rather than the bridge's numeric native gene identifiers.

Checks include the exact candidate universe and original candidate fields,
parent identity and duplication flag, parent degree, sister gene/taxon counts,
modeled nonfocal coverage, eligibility priority, every nearest tie, lexical
representative and all four sequence distances. Tie membership must match
exactly under the original absolute 1e-12 rule. Numeric distances allow 1e-12
absolute/relative rounding differences between downward summation and upward
`math.fsum`; maximum observed differences are reported. This validates the
selection algorithm on fixed trees, not orthology, rooting or event timing.

The full run completed successfully under
`metadata/duplication_sister_reference_readback_plan_20260926.json`, with its
process identity in the corresponding launch record. All inherited source
pins were checked equal to the original inventory plan. Resources are one
CPU, 8 GiB RAM, no swap and negligible output. The 0.1–8 hour planning range
is uncalibrated. Run with:

```bash
python scripts/readback_duplication_sister_references.py \
  --plan metadata/duplication_sister_reference_readback_plan_20260926.json
```

Output is `results/orthology/duplication-sister-reference-readback-20260926-v1/`.
The final receipt reports `passed_full_duplication_sister_reference_readback`.
All 109,245 profile and 109,228 MAFFT rows passed, with exact reference
choices, ties, counts and eligibility statuses. Maximum sequence-distance
difference was 3.552713678800501e-15 in each guide. All source pins and output
artifact hashes were rechecked after successful service exit. The receipt is
archived in `metadata/duplication_sister_reference_readback_completed_20260926.json`.
Known-tree tests passed all six statuses, ties, four known distances,
isolation from a very long ancestral edge, six deliberately corrupted export
fields and malformed-tree rejection. Reproduce these with
`python scripts/check_duplication_sister_reference_readback.py`.

## Additional alignment-input preparation queued

`scripts/materialize_duplication_reference_inputs.py` is queued behind the
exact supplementary coordinate-readback process. It requires a successful
readback bound to the coordinate producer and the reference inventory, checks
the exact additional-model partition again, and verifies every source record's
model, version, raw-file path/hash and sequence hash. Every coordinate shard
must have a matching independent proof. It produces full and pLDDT>=70
C-alpha PDBs using the same serializer as the primary workflow, preserving
original residue positions, masked sequences and file hashes. Short masks and
rejected source models remain explicit dispositions.

The planned scope is 14,540 models and 29,080 model/mask dispositions. Resources
are one CPU, 8 GiB RAM, no swap, a 4 GiB output planning allowance and a 100 GiB
free-disk reserve; uncalibrated runtime range 0.1–8 hours after the readback.
The job uses no GPU or paid services. Configuration and exact process identity
are in `metadata/duplication_reference_alignment_input_plan_20260926.json`
and the corresponding launch record. Reproduce with:

```bash
python scripts/materialize_duplication_reference_inputs.py \
  --plan metadata/duplication_reference_alignment_input_plan_20260926.json
```

Output will be
`results/structural_comparisons/duplication-reference-alignment-inputs-20260926-v1/`.
Its completion status explicitly identifies additional reference inputs.
A synthetic completed handoff passed full, sparse-mask, empty-mask and
source-rejected cases, written-coordinate hashes and altered-proof rejection:
`python scripts/check_duplication_reference_materialization.py`.

Production input preparation is waiting, not complete. Reference-pair
alignment orchestration still needs to combine these additional inputs with
the relevant audited primary inputs. Full and masked scores normalize to
different retained lengths; no PAE or domain-orientation qualification is
implied by serialization.

## Complete additional reference alignment run queued

`scripts/run_duplication_reference_alignments.py` now waits for both exact
input-preparation processes. The handoff helper checks successful plan-bound
receipts, exact model/mask universes for both collections, source hashes,
manifest hashes and all disposition counts. It rechecks the reference
inventory's independent ledger binding and exact additional-model partition.
Canonical pair hashes and membership in the frozen primary pair table decide
which pairs require additional computation. All 321 already-covered primary
pairs remain in the inventory for later result joining; they are not treated
as completed here or rerun by the supplementary job.

The full workload is 32,541 additional pairs × two input orders × two masks,
or 130,164 explicit dispositions. The job reuses the unchanged primary
per-pair runner, preserving raw native outputs, timings, command/input hashes,
checkpoint bindings and explicit unavailable/error/timeout outcomes. A digest
of both materialization manifests and receipts binds every checkpoint.

Resources are two CPU workers, 8 GiB RAM, no swap, a 24 GiB output planning
allowance, 100 GiB free-disk reserve and 600-second timeout per native call.
The 8–480 hour planning interval is scaled from the primary uncalibrated
range, not a measured ETA or guaranteed bound. No GPUs or paid services are
used. At most 64 jobs are pending in the worker pool. Configuration and exact
process identity are recorded in
`metadata/duplication_reference_alignment_plan_20260926.json` and the
corresponding launch record. Reproduce with:

```bash
python scripts/run_duplication_reference_alignments.py \
  --plan metadata/duplication_reference_alignment_plan_20260926.json
```

Output will be
`results/structural_comparisons/duplication-reference-alignments-20260926-v1/`.
The successful completion status will explicitly retain `pending_readback`:
this is not independent numeric verification or evidence of biological
asymmetry. Joining reused primary pairs, numeric verification, confidence and
coverage interpretation, domain/orientation controls and statistical tests
remain pending.

`python scripts/check_duplication_reference_alignment_handoff.py` passed a
native two-source completed-handoff fixture: the primary pair was excluded,
the additional pair aligned in both directions, short-mask outcomes were
preserved in both directions, all checkpoint hashes matched, and altered
source provenance was rejected even after updating the manifest hash.
Production alignments are queued, not yet started.

## Additional coordinates verified and alignment inputs prepared

The supplementary raw-coordinate validation and independent readback completed
for all 14,540 models, with no rejections. Every exported accepted coordinate,
sequence and confidence value was reconstructed from the raw CIF atom rows:
6,029,194 C-alpha residues in total. All source pins, shard output hashes and
independent proof hashes were checked again. The result is archived in
`metadata/duplication_reference_coordinates_completed_20260926.json`.

Full/masked input preparation also completed with 29,080 dispositions:
14,540 full inputs ready, 14,050 pLDDT>=70 inputs ready, and 490 masked inputs
with fewer than three retained residues. These are explicit exclusions, not
failed predictions or missing structures. The 28,590 ready PDBs total
816,376,377 bytes. The exact model/mask grid, all ready PDB hashes, byte totals,
source pins and upstream receipt bindings were checked and archived in
`metadata/duplication_reference_alignment_inputs_completed_20260926.json`.
That check is a manifest/hash readback, not an independent numeric validation
of every serialized PDB; the queued alignment readback handles the latter for
successful comparisons. Reference alignments still wait for primary inputs.


## Sequence covariates for the asymmetry analysis

The complete provisionally referenced cohort now has fixed sequence-tree
covariates in `results/orthology/duplication-sequence-covariates-20260926-v1/`.
These preserve every source field and add terminal sequence lengths for each
copy, their signed difference, and the difference divided by duplicate-pair
sequence distance. Gene A/B orientation is lexical and fixed between guides;
its sign has no intrinsic biological meaning. Tip lengths are obtained by
subtracting the common duplicate-node-to-reference path from each tip-to-reference
path. The pair sum is checked against the independently exported duplicate path.

Numerical tolerance is 1e-12 times the larger of one and the source path lengths.
Negative residuals only within tolerance are clamped and counted; incompatible
paths are rejected. Pair distances at most four tolerances have blank normalized
contrasts, and directional differences at most two tolerances are unresolved.
This is a floating-point resolution rule, not an uncertainty interval or a
biological effect-size threshold. Sequence lengths are relative divergence,
not change per calendar time. Reference eligibility and gene-tree uncertainty
remain as documented above.

| Sequence covariate status | Profile | MAFFT |
| --- | ---: | ---: |
| Provisionally referenced events | 17,461 | 17,448 |
| Resolved pair-distance denominator | 15,868 | 15,858 |
| Near-zero pair-distance denominator | 1,593 | 1,590 |
| Unresolved direction (including near-zero pairs) | 1,750 | 1,747 |

Of 17,392 shared pairs, 15,643 have the same resolved direction across guides,
six have opposite resolved directions, and 1,743 are unresolved in one or both.
The full comparison also retains reference-set agreement. These are covariates
for future structure comparisons, not independent replicates, structural
asymmetry results, significance tests or evidence of positive selection.

All 34,909 output records retained exact original source fields. Independent
decimal three-tip algebra reconstructed terminal lengths and contrasts, and
all 17,392 cross-guide comparisons were checked. Artifact hashes and the six
direction-changing cases are archived in
`metadata/duplication_sequence_covariates_completed_20260926.json`.
Known-distance tests cover sign reversal, scaling, equal/zero/near-zero paths
and incompatible/nonfinite inputs. Reproduce with:

```bash
python scripts/check_duplication_sequence_covariates.py
python scripts/prepare_duplication_sequence_covariates.py \
  --inventory results/orthology/duplication-sister-reference-inventory-20260926-v1 \
  --readback results/orthology/duplication-sister-reference-readback-20260926-v1 \
  --output <fresh-output-directory>
```

This small table transformation uses one CPU and no GPU, with no new native
phylogenetic fitting. Structural-response joins, sequence-conditioned tests,
family/phylogenetic dependence and nonduplication controls remain unfinished.

## Primary whole-protein numerical qualification (September 27)

Full primary geometry and input-order readbacks now pass. All 387,646 successful
native mappings were reconstructed and checked; 327 short mappings are
numerically degenerate, exactly matching the separate one-/two-residue census.
The original strict audit failure remains preserved. All 27 native RMSD rounding
discrepancies occur among these excluded two-residue mappings.

The complete order-summary readback verifies 206,400 pair/mask rows and 4,647,688
numeric values, retaining all 103,200 primary model pairs:

| Numerical disposition | Full protein | pLDDT ≥70 masked input |
|---|---:|---:|
| Both input orders usable | 103,200 | 90,442 |
| One input order usable | 0 | 35 |
| Both orders numerically excluded | 0 | 146 |
| Both orders lack usable inputs | 0 | 12,577 |

The 35 asymmetric cases comprise nine with only order 0 usable and 26 with only
order 1 usable. Neither order is selected to rescue the pair for two-order
sensitivity analysis. The 146 double exclusions plus 35 single exclusions
account for all 327 excluded directions. Missing masked inputs account for
25,154 directed dispositions. Thus masking changes measurement availability as
well as the coordinates being compared; pair membership must be held constant
for a direct mask-sensitivity contrast.

These are numerical eligibility counts, not independent evolutionary events or
coverage-qualified biological comparisons. Source controls, aligned fraction,
domain orientation, prediction uncertainty, phylogenetic dependence and matched
backgrounds remain relevant. The residue-correspondence sensitivity stage now
examines all 193,642 pair/mask comparisons having two usable input orders, with
90,442 pairs shared across the two masks. Its independent check remains pending.

Source-bound closure records:
`metadata/primary_diagnostic_geometry_completed_20260927.json` and
`metadata/primary_order_summary_completed_20260927.json`.

## Whole-protein common-residue stage launched

The full 17,619 oriented A/B/reference model triads now expand to 281,904 records
(two masks × eight independent edge-order choices), preserving all 36,944
associated event/reference links. For each combination, the mapper stores the
original-position triples shared through the reference, and separately the
subset consistent with the direct A–B correspondence. This prevents different
pairwise residue subsets from being silently treated as a common structural core.
All missing and numerically excluded edge states remain explicit.

There are 15,713 three-model, 1,904 two-model and two one-model triads. Repeated
model identities use explicit identity correspondences; they are not new native
alignments or independent structural observations. Order combinations remain
explicit even when an identity edge makes some computations repeat. Coordinate
fits and their geometry/coverage qualification are later stages.

Producer and independent full readback use
`metadata/whole_protein_common_residues_plan_20260927.json`. Full serialized
verification remains pending; no asymmetry result is claimed.

The common-coordinate fitting and independent quaternion readback are queued
behind the full map audit. Both reference-common and cycle-consistent definitions
are retained (563,808 fit rows). Coverage uses each original full-protein length;
30/50-residue minima crossed with 0.5/0.7/0.9 coverage preserve the existing six
screen settings. These screening outputs remain conditional measurements, with
shared-model identity and all source exclusions explicit. Signed differences
compare distances on identical residue triples; they are not reconstructed
ancestral changes or tests of asymmetric evolution.

## Complete whole-triad mapping independently verified

Both mapping services terminated successfully and all 281,904 records passed
independent reconstruction. They contain 58,204,472 reference-common residue
occurrences and 53,654,022 cycle-consistent occurrences. These repeat across
masks, orders and linked events; they are not independent residue observations.

There are nonempty reference intersections in 262,912 records and nonempty
cycle-consistent subsets in 228,566. The two definitions differ in 175,720
records, underscoring the need to retain both rather than assume three pairwise
alignments are transitive. Missing and numerically excluded edges remain explicit.
Full coordinate fitting has started after successful readback; its geometry,
coverage screens and independent numerical verification remain pending.

Closure evidence: `metadata/whole_protein_common_mapping_completed_20260927.json`;
full checker output: `metadata/whole_protein_common_residues_readback_20260927.json`.

## Expanded catalog reference pipeline completed September 30

Full expanded reference choices and their independent native-tree/path
reconstruction have completed for **283,409 modeled duplicate target links**
across 50,063 candidate-bearing trees. Every unavailable or ineligible target
remains explicit:

| Native guide | Targets | Provisional reference | No modeled nonfocal sister | Parent reported duplication | Nonbifurcating parent |
|---|---:|---:|---:|---:|---:|
| Profile | 141,724 | 29,155 | 40,700 | 71,765 | 104 |
| MAFFT | 141,685 | 29,130 | 40,695 | 71,757 | 103 |

The independent check uses native leaf-interval subtraction and upward
`math.fsum` tree paths; maximum path differences are 3.55e-15 under each guide.
All nearest modeled tied references remain in the ledger, without selecting
on structural responses. The full **121,490 event/reference/duplicate-side
rows** involve 83,207 models and 55,701 distinct model pairs, including 55,336
pairs additional to the expanded primary queue. This is a queue partition,
not a count of unmeasured pairs after earlier reference/background reuse.
The additional model partition has 24,804 models and 9,597,401,251 stat-counted
raw coordinate bytes; full coordinate parsing is a separate requirement.

A stronger independent reader verifies every gene against the frozen native
bridge, every model/version/sequence/path assignment, full exported catalog
records, all tied-reference/side combinations and pair/work partitions.
It checked **86,440 native gene/model mappings**. Isolated fixtures include
numeric versions 6 and 10 and reject rehashed but incorrect gene/model
assignments, full catalog metadata changes, false work partitions and missing
tied sides. Original producers and earlier readers remain unchanged.

Independent cross-guide outer-join readback covers all **142,107 gene pairs**:
141,302 shared, 422 profile-only and 383 MAFFT-only. Among 29,016 provisionally
eligible under both guides, **28,950 have identical nearest-reference gene
sets** and **66 have disjoint sets**; none partly overlap. Chosen model
identity agrees for 28,951, one more than chosen gene identity, because distinct
reference genes can share a model. Model agreement is not gene identity or
independent prediction agreement.

[Completion evidence](../metadata/expanded_duplication_reference_pipeline_completed_20260930.json)
binds 82 source/artifact hashes and all six exact captured-process completion
journals. These are extant, availability-dependent provisional references,
not ancestral states or demonstrated biological orthologs. Native reference
orthology, sequence-locked choice sensitivity, coordinate/input reuse and
quality, confidence/PAE/domain controls, common-residue triads, phylogenetic and
ancestral uncertainty, and calibrated asymmetry tests remain required. Earlier
common-triad outputs retain their older frozen scope. None of the eight
scientific aims is complete.

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

Full expanded reference input preparation and independent written-PDB checking
are queued after successful raw-coordinate readback for all 24,804 additional
models and both masks (49,608 dispositions). The v2 materializer preserves the
serialization algorithm while requiring the stronger native-ledger partition
loader. Every short/rejected input remains explicit; shared models across
background/reference inputs require exact mask/position/byte equality before
union reuse. Full end-to-end fixtures reject rehashed wrong coordinates,
versions and residue mappings. No supplemental production completion, native
alignment or structural asymmetry result is claimed. Plans:
`metadata/expanded_duplication_reference_alignment_input_plan_20260930.json`
and `metadata/expanded_duplication_reference_alignment_input_readback_plan_20260930.json`.

## Sequence-first reference choice completed September 30

The complete sensitivity analysis selected nearest nonfocal sister genes from
sequence-tree path distances **before** looking up structures. It retained
all **283,409 modeled duplicate targets** across **50,063 native family trees**,
every original available-reference field, all nearest ties within 1e-12, a
fixed lexical representative, and missing model/version assignments. Choosing
a modeled member of a tie does not silently replace an unmodeled lexical
representative. Every parent-ineligible target stays explicit and is not
promoted to an acceptable biological reference.

Independent native leaf-interval subtraction and upward `math.fsum` paths
checked every saved field and the full original-context/choice matrices.
Maximum path disagreement is **3.55e-15** under each guide. The completion
record binds **59 source/artifact hashes and two exact captured-process
journals**:
[full evidence](../metadata/expanded_sequence_first_references_completed_20260930.json).

| Sequence-first disposition | Profile targets | MAFFT targets |
|---|---:|---:|
| Closest lexical gene has a model | 21,581 | 21,549 |
| Closest ties have no models; available-reference design uses a farther gene | 7,560 | 7,566 |
| Lexical gene unmodeled but another equally nearest tie has a model | 14 | 15 |
| No modeled nonfocal sister reference | 40,700 | 40,695 |
| Original parent context ineligible | 71,869 | 71,860 |
| All modeled duplicate targets | 141,724 | 141,685 |

Among the **29,155/29,130 provisionally available reference targets**, the
farther-gene category is **25.93%/25.97%**. These are conditional record
fractions, not independent event rates, uncertainty intervals or evidence of
structural asymmetry. The two guide counts overlap. This result demonstrates
that structure availability changes reference distance in this frozen design;
it does not establish whether that changes a structural contrast. The full
per-target distance increments and original/sequence-first model choices are
retained, including ineligible contexts and unavailable genes.

This sensitivity complements the existing available-reference ledger; it
neither replaces it nor rematches fixed background controls. Extant nearest
genes are not ancestors, and a parent not reported as duplication is not proof
of speciation. Reference orthology and gene-tree/root/model uncertainty remain
required. Coverage/confidence/PAE/domain controls, common-residue mapping, structural
measurement reuse and
calibrated effect tests are still unfinished. Missing sequence-first models
remain missing until suitable source-qualified structures become available;
no predictor is silently switched and no GPU inference has launched.

Known-case and full-handoff fixtures cover nearer unmodeled genes, model
versions 6/10, modeled/unmodeled ties, no modeled references, parent/focal/
multifurcation exclusions and isolation from very large outer-ancestor branches.
Four rehashed false choice/distance/model/status exports are rejected.
Both production and independent reader used one CPU/16 GiB/no swap. Resource
planning predates launch; 0.1–8-hour ranges were uncalibrated, not ETAs. Scripts,
plans, input/output provenance and exact launches are versioned in
`metadata/expanded_sequence_first_references_queued_20260930.json`.
This gene-choice check is distinct from the older sequence-locked *residue*
correspondence control. All eight scientific aims remain incomplete.

## Expanded reference coordinates and written inputs independently verified

Full native-CIF readback finished for all **24,804 additional reference models**,
**9,982,283 residues** and **25 shards**, with no content rejections. It used
9,597,401,251 raw coordinate bytes and preserved the exact model partition.
The 77-source/artifact completion closure includes every coordinate shard and
its full independent proof plus both exact captured-process journals:
[coordinate evidence](../metadata/expanded_duplication_reference_coordinates_completed_20260930.json).

Full supplemental PDB input preparation and independent checking also finished
for **49,608 two-mask dispositions**, including **24,804 full-protein ready
inputs**, **23,781 pLDDT70 ready inputs** and **1,023 confidence masks with fewer
than three retained residues**. Short masks remain excluded inputs, not zero
structural distances. Every original residue position, amino acid, model/
version, source checksum, written coordinate/confidence rounding, PDB hash and
complete source/shard/proof lineage was checked. Verified PDB bytes total
**1,348,146,870**; completion evidence is
`metadata/expanded_duplication_reference_inputs_completed_20260930.json`.

These completed inputs advance the full 55,701 distinct-pair reference ledger;
they do not establish complete pair measurements or authorize reuse of old
results. The prior catalog-source screen still requires raw/materialized input,
mask/position, executable/settings, input-order, checkpoint and numerical proof
checks before reuse. Full numerical/coverage/confidence/PAE/orientation quality,
common-residue geometry, reference orthology and biological asymmetry remain
pending. Background input preparation still waits for its complete coordinate
readback. GPU inference remains paused; all scientific aims remain incomplete.


## Full sequence-first guide comparison completed September 30

Independent dataframe outer merging checked every exported field, matrix and
summary for the complete **142,107 duplicate gene-pair union**: **141,302 shared
pairs**, **422 profile-only** and **383 MAFFT-only**. Original gene order is
canonicalized only for identity; each guide's family, context, nearest ties,
lexical gene, missing model/version and distance fields remain separate.
All 283,409 native target rows and their completed sequence-tree proof remain
bound. Full closure checks **76 source/artifact hashes and two exact process
completion journals**:
[completion evidence](../metadata/sequence_first_guide_comparison_completed_20260930.json).

| Comparison within the full pair union | Pairs |
|---|---:|
| Eligible sequence-sister context under both guides | 69,683 |
| Same nearest sequence-gene set and lexical gene | 69,535 |
| Disjoint nearest sequence-gene sets | 148 |
| Overlapping but unequal nearest sets | 0 |
| Parent context ineligible in one or both guides | 71,619 |
| One-guide candidates | 805 |
| Both lexical references have frozen models, among eligible shared targets | 21,451 |
| Same chosen model/version, among those modeled targets | 21,411 |

Sequence eligibility includes eligible parents with no modeled sister, so this
comparison is not restricted to the older available-reference subset. Two
missing model assignments never count as agreement between predictions; model
agreement is evaluated only when both choices have models. Different reference
genes can share a model, and model agreement does not establish gene orthology.
No missing lexical gene is silently replaced by an equally near modeled tie.
The paired guide records overlap and do not represent independent events.

Nine known cases and a complete synthetic two-guide handoff passed, including
reversed duplicate labels, disjoint/overlapping ties, modeled and unmodeled
choices, shared models from different genes, parent exclusions and one-guide
candidates. Five rehashed false exports were rejected. Producer and independent
reader each used one CPU/8 GiB/no swap; no native alignments or predictions.
This checks guide-choice sensitivity, not accepted phylogeny, reference
orthology, ancestral states or structural asymmetry. Sequence-first structural
contrasts and missing-reference recovery remain outstanding.

## Full new reference measurements launched September 30

The audited inputs now feed native USalign measurements for all **24,947 current
reference pairs without matching existing catalog sources**, at most **99,788
mask/order dispositions**. Both input orders use the original full-protein and
pLDDT70 masks, executable and options. Every successful native output and every
short/error/timeout disposition is checkpointed; failed calls are not silently
retried or replaced. This is the complete new-work partition of the full design,
not a pilot or a selection based on structural responses.

The full **55,701-pair work ledger** also retains **30,754 pairs with matching
prior catalog sources** as pending actual input/result/numerical reuse checks.
No old result is imported or numerical flag cleared by this partition. Catalog
checksums alone do not prove compatible masks, residue positions, native input
bytes/settings, checkpoints or numerical quality. The final result union must
cover all 55,701 pairs and preserve all earlier exclusions. Seventeen new
reference pairs also occur in the expanded background workload; future union
assembly must check their exact inputs/results before sharing measurements.

A full actual-data handoff preflight checked **602,972 two-mask input
dispositions**, exact current model partition, written-input proofs and all
reference/source partitions before launch. A complete synthetic native handoff
and diagnostic passed both input orders, short masks and rehashed false
checkpoint, incomplete-partition and wrong-source rejection. Existing scripts
and datasets remain intact; versioned adapters handle the stronger native
reference audit and the full source partition.

Resources were estimated before launch: **eight single-thread native workers,
32 GiB RAM, no swap, 24 GiB output allowance and a 100 GiB free-disk reserve**.
A deterministic timing sample of the completed reference run has a mean native
elapsed time of 0.474 seconds; a simple eight-worker throughput projection is
1.64 hours. New endpoint lengths span 29–1,280 residues (median 319). The
**2–48-hour planning range is not an ETA**: differing lengths, scheduling,
batching, failures and tails limit extrapolation. The sample uses existing
measurements for resource planning and does not run a new protein pilot.
GPU predictions remain paused; no paid resources are used.

Exact dependency wrappers also queue the full native-text/residue/least-squares
RMSD diagnostic (retaining all discrepancies), coordinate-rank/rotation-curvature
assessment and independent quaternion geometry reader. They use two CPU/32 GiB
for the diagnostic and one CPU/16 GiB per geometry stage, no swap. These stages
are queued, not complete. Plans, resources and exact process identities:
`metadata/expanded_reference_alignment_pipeline_queued_20260930.json`.

Full result reuse/union, input-order correspondence, original-length coverage,
common-residue triads, sequence-locked residue controls, confidence/PAE/domain
orientation checks, reference orthology and phylogenetic/calibrated asymmetry
inference remain required. None of the eight scientific aims is complete.


## Actual reference input/result reuse qualification — September 30

The full reference reuse gate now has independently verified outputs for all **55,701 distinct
pairs**, **two masks and two input orders**, totaling **222,804 dispositions**.
It checks every eligible old source against the current frozen inputs rather
than treating catalog checksum matches as permission to reuse results. All
**24,947 new pairs / 99,788 dispositions** remain explicitly pending the
currently running native measurements; they are not filled from old results.

For the **30,754 matching-source pairs**, a deterministic source preference
uses the completed expanded primary run for **365 pairs** and the earlier
reference run for **30,389 pairs**. This preference uses the completed source
role, not structural responses or numerical usability. Other matching catalog
source labels remain in the original source-union ledger. A better-looking
alternative run cannot silently clear the selected source's numerical flags.

Qualification compares raw coordinate/sequence checksums and original length;
actual raw bytes; mask and disposition; every original residue position and
retained amino acid; actual old/current PDB bytes; native executable checksum
and options; complete directed checkpoint provenance and native text; and the
full numerical and rigid-fit geometry records with their independent original
proofs. Source collection paths/shard numbers may differ, but coordinate
content and residue mappings must agree. Original checkpoints remain unchanged,
with their paths and checksums retained. Original input order is mapped by
actual ordered model/version endpoints, including reversed source directions.
No zero distances replace short, missing or failed measurements.

The producer found **123,016 compatible old directed dispositions**, with no
incompatible input signatures. It retained **118,804 numerically usable
alignments**, **4,154 unavailable inputs**, **53 short/nonunique rotations**,
**four short/nonunique mappings also flagged for RMSD discrepancy**, and
**one additional RMSD discrepancy**. These counts passed the full independent reader and exact-process completion
closure. They count directions,
not independent biological events, unique pairs or accepted asymmetry effects.
A compatible old result is not automatically a usable numerical measurement.

Checks cover **46,200 unique reused endpoint models**, **17,405,911,378 raw
coordinate bytes by initial stat**, and **92,474 model/mask/source checks**;
the latter retain models used by both preferred source roles rather than
assuming independent observations. The producer's complete byte/provenance
index has **348,736 source bindings**. Large exported state, proof and byte
indexes remain under
`results/structural_comparisons/reference-reuse-qualification-20260930-v1/`,
outside Git. Versioned plans/scripts record their locations and source pins.

The independent reader reconstructs every exported input identity and all
222,804 directed states from original source manifests/checkpoints and complete
numeric/geometry proof tables. It separately derives source preference,
canonical input signatures, directed correspondence and numerical exclusions;
source I/O/proof loading and native-text parsing libraries are shared. It
rehashes raw/PDB/checkpoint bytes and checks every total. This is export/result
qualification, not a new native optimization or an independent predictor.
The separate exact-dependency closure also finished: **348,748 source/artifact
hashes and both captured process completion journals** passed. Its full byte/
checkpoint index remains outside Git, with the versioned location/checksum in
[the completion locator](../metadata/reference_alignment_reuse_completed_20260930.json).
The independent reader and completion closure are complete; no old numerical
exclusion was cleared.

A complete native synthetic handoff passed all **16 states**, including model
versions 6/10, reversed old directions, both masks, new pairs, incompatible
masks, unavailable inputs and nonunique rotations. Six rehashed false exports
were rejected, including cleared exclusions and false native/source/input
identities. Synthetic proof stubs are test data only; they do not qualify any
production protein. This is a full-workflow check, not a protein pilot.

Resources were recorded before launch: two CPU/16 GiB/no swap for the producer
and independent reader, eight GiB output allowance each, an uncalibrated
0.1–8-hour planning range, and no GPU or paid resources. Closure uses one
CPU/16 GiB/no swap. Exact plans and handles are in
`metadata/reference_alignment_reuse_queued_20260930.json` and
`metadata/reference_alignment_reuse_closure_queued_20260930.json`.

The complete independent reuse proof is now available; the full union with
new measurements is still required. Input-order correspondence, original-length coverage,
common-residue triads, sequence-locked correspondence, domain/orientation/PAE
and prediction controls, native reference orthology, phylogenetic uncertainty
and calibrated biological asymmetry remain unfinished. All eight scientific
aims remain incomplete.

## Full measurement union queued — September 30

The complete reference design has **55,701 pairs**, two masks and two input
orders, totaling **222,804 dispositions**. The union workflow preserves all
**123,016 independently qualified old states** and waits for all **99,788 new
states** rather than exporting an apparently complete table with pending
measurements silently omitted. Old sources retain the existing preference
for expanded primary results, otherwise old reference results; numerical
usability never selects an alternative source or changes a flag.

Four additional stages are now launched with exact captured dependencies:

1. Complete new-native measurement verification binds every checkpoint and
   checks the complete pair/mask/order grid, numerical and quaternion proof
   lineage, every disposition and all four original process-completion
   journals. Its large checkpoint index stays outside Git.
2. The union exports every current pair/mask/order with actual ordered
   model/version endpoints, original checkpoint path/hash, native source
   order, native metrics and complete original numerical/geometry records.
   Old native source order can differ from current order; ordered endpoints
   establish the correspondence. Original checkpoints are unchanged.
3. A separate reader independently reconstructs every union field and all
   source/status/exclusion counts from the completed source ledgers and
   actual checkpoint bytes. It shares source I/O/proof loading, not producer
   projection or exclusion logic.
4. Final closure requires the full independent proof, all exported/raw source
   bindings and both union producer/reader completion journals. A small
   versioned locator records the outside-Git full proof archive/checksum.

No missing, short, timeout, parse/native error or degenerate/RMSD-flagged result
is replaced by zero. No native optimization or rigid-fit measurement is
repeated in this union. Numerical usability is a provenance/identifiability
label and remains distinct from original-length coverage, confidence,
common-residue/sequence-locked correspondence, domain orientation, PAE,
prediction-source and biological orthology/phylogenetic/calibration controls.
Masked native coverage is not silently relabeled as full-protein coverage.

The synthetic handoff passed all **16 states** across four pairs, including
model versions 6/10, reversed old-primary directions, both masks, unavailable
inputs, nonunique geometry and explicit synthetic RMSD/error/parse/timeout
states. Six rehashed false exports were rejected: cleared exclusions, wrong
source order, reversed endpoints, changed native metric, wrong checkpoint
hash and deleted unavailable state. Source/numeric/quaternion/journal proof
stubs are fixture data, not production qualification or a pilot. The actual
completed reuse proof schema, full scope, numerical counts and target native
plan were also checked before queueing.

Resources were recorded before launch: two CPU, 16 GiB memory, no swap,
eight GiB output allowance and an uncalibrated 0.1–12-hour planning range per
stage, with a 100 GiB disk reserve checked before launch. These stages wait
behind the existing native pipeline; no GPU or paid resource is used. Exact
plans and live handles are in
[the queue record](../metadata/reference_measurement_union_pipeline_queued_20260930.json).
The existing native measurements remain active. These stages are queued;
their production outputs and completion proofs are not yet present.

The workflow check is reproducible with
`python scripts/check_reference_measurement_union.py`. Production plans pin
all scripts and source plans; their output directories are append-only and
completed data are never overwritten. Full pair/order summaries, coverage
and event projection, common-residue triads, sequence-locked/domain/PAE/
prediction controls and calibrated phylogenetic asymmetry remain downstream.
All eight scientific aims remain incomplete.

## Full native reference orthology query — September 30

The new native-membership workflow uses all **283,409** independently checked
target contexts, retaining every original field and both the availability and
sequence-first reference designs. It includes all nearest ties, unmodeled
reference genes, lexical choices, absent references and parent-ineligible
contexts. The full preflight counted **214,461 reference tie records** and
**428,922 duplicate-to-reference links**. These overlap across designs and
guides and are not independent events. This universe is larger than the
121,490-row modeled, parent-eligible structural comparison ledger.

Profile-guide availability/sequence-first ties total **32,533/74,720**;
MAFFT-guide ties total **32,507/74,701**. Original parent-ineligible contexts
contain **2,148/2,147 availability ties** and **2,751/2,759 sequence-first ties**.
Their native memberships are recorded for sensitivity checks; a positive
membership never promotes an excluded parent. Empty gene sets remain explicit
context records and are not coded as absent native orthology.

Each reference is queried against both duplicate genes in both audited native
orthology outputs. Identical physical gene-pair queries are deduplicated while
every context/design/tie/side link remains. The compiled merge scans both
complete reciprocal streams (**31,448,413,556 bytes** combined), checking
global ordering, canonical keys and exactly one record per direction.
Source-line protein ordinals are bound through the grouped-protein identity
audits and the completed resolved-tree membership proofs. The profile v3 path
resolves to the original v2 profile directory; both source mappings also have
identical checksums. No sequence-tree family/version substitution is assumed.

The separate reader reconstructs the whole context/link ledger, model/version
assignments, nearest ties and parent exclusions. It checks every unique native
query by fixed-record binary search, with reciprocal validation and bracketing
for absence. It separately derives whether a reference is assigned to both
duplicates, only one, or neither under each native guide. Source I/O/proof
loading and the fixed-record reader are shared; producer ledger and merge
logic are not. Full byte bindings and both exact completion journals are
required before completion is reported.

Synthetic adapter checks passed the complete ten-context fixture and rejected
six rehashed false exports: cleared parent exclusion, changed lexical choice,
deleted unmodeled context, changed cross-guide membership, deleted nonlexical
tie and changed native ordinal. The fixture uses stub native-tree/source
qualification proofs and does not qualify a production protein. The initial
adapter failed before processing because its expected source status omitted
the word `full`. Failed scripts/plans/journals are preserved. Fresh v2 scripts
and outputs correct only that label; all biological criteria and source data
remain unchanged. Actual production preflight then passed **86 source
bindings**, including complete stream hashes and native ordinal lineage.

Resources were recorded before launch: one CPU, 16 GiB memory, no swap, four
GiB output allowance and an uncalibrated 0.1–8-hour planning range for each
processing/checking stage. Closure uses one CPU/16 GiB/no swap. No GPU or paid
resource is used. Exact plans and live launch handles are in
[the v2 pipeline record](../metadata/reference_orthology_queued_20260930_v2.json).
Native membership, independent readback and journal closure are complete:
**166,829 distinct gene-pair queries**, all **428,922 logical links**, **103
source/artifact bindings** and **two exact process completion journals** passed.
The profile/MAFFT streams contain **157,740/157,754** queried pairs and lack
**9,089/9,075** respectively. These totals include parent-ineligible contexts
and both reference designs; absence is native output absence, not a validated
biological loss. See [the full completion proof](../metadata/reference_orthology_completed_20260930_v2.json).

The full premeasurement cohort summary also passed independent SQL
reconstruction of all **566,818 context/design records**, **48 nested policy
rows** and **76 lexical reference/context/missingness/native-state cells**.
Its closure checked **23 bindings and two additional exact completion
journals**. The complete [policy table](tables/reference_native_orthology_context_policy_counts_20260930.tsv)
and [context-state matrix](tables/reference_native_orthology_lexical_reference_context_matrix_20260930.tsv)
are exact copies of the generated, independently checked source tables.

| Native source guide | Reference design | Parent-eligible contexts | Lexical reference with a model | Lexical model assigned to both duplicates in both native guides |
| --- | --- | ---: | ---: | ---: |
| Profile | Availability | 69,855 | 29,155 | 28,569 |
| MAFFT | Availability | 69,825 | 29,130 | 28,551 |
| Profile | Sequence-first | 69,855 | 21,581 | 21,169 |
| MAFFT | Sequence-first | 69,825 | 21,549 | 21,150 |

Among parent-eligible sequence-first contexts, **68,835/68,829** lexical
references are assigned to both duplicates under both native guides. Of
those, **47,666/47,679** lack a model in the frozen catalog. Missing structures
are therefore recorded separately from native orthology disagreement. These
are conditional source-context counts, not independent events or final
structural sample sizes. Any/all tied-reference native policies agree in this
catalog; that observation never authorizes replacing a lexical unmodeled
reference with another modeled gene or a farther reference. Biological
orthology and structural qualification remain separate requirements.

These checks establish native output membership, not independent biological
orthology, accepted tree rooting, ancestral structure or a structural
duplication effect. Small-family orthology supplements remain excluded from
the native streams. Full structural measurement union, sequence-first
structural sensitivity, original-length coverage, common-residue and
sequence-locked comparisons, domain/PAE/prediction controls, phylogenetic
uncertainty and calibrated inference remain required. All eight scientific
aims remain incomplete.

## Complete context-to-measurement design — September 30

The full native source now has a verified fixed model/version work assignment
for all 283,409 target contexts under both reference designs. The complete
214,461 tied-reference records/428,922 logical duplicate-reference sides are
retained, including missing reference structures, parent-ineligible contexts
and identical models. All 121,490 original availability-side ledger links
were checked. Duplicate gene identity determines each frozen model/version
role even when queue and reference table positions differ. Model pair hashes
and current focal endpoints are explicit; reference gene aliases sharing a
model remain separate logical links.

Completion evidence is
`metadata/reference_context_measurement_design_completed_20260930.json`: all
40 source/artifact hashes and both exact original process journals checked.
`metadata/reference_context_measurement_design_completed_readback_20260930.json`
independently reconstructs every original source/native field, duplicate gene/
model/version mapping, reference design/tie/lexical choice, parent exclusion,
physical pair endpoint and availability link using uniquely keyed SQL joins.
Both complete context and availability-link universes and all native source
ordinals were exhausted. The large immutable export remains outside Git in
`results/orthology/reference-context-measurement-design-20260930-v1`.

[Verified work-disposition counts](tables/reference_context_measurement_work_dispositions_20260930.tsv)
label counting units explicitly. Empty-reference cells count contexts; other
cells count duplicate/reference side links. Parent-eligible distinct sides in
the physical full-reference design number 60,406/60,354 for availability and
44,668/44,608 for sequence-first designs (profile/MAFFT). These overlap across
guides, designs and ties and do not count unique predictions, independent
events or accepted structural cohorts. Missing sequence-first structures remain
explicit, including an unmodeled lexical gene with a modeled tied alternative.

One parent-excluded profile side is also in the physical measurement design
under each reference design. Those two design views overlap and remain
excluded. Presence of a structural comparison for another event never
overrides the original parent context. Subsequent full-context coverage and
own-guide/both-guide/lexical/any/all-tie policies must enforce that exclusion.

Reproduction uses `prepare_reference_context_measurement_design.py`, then
`readback_reference_context_measurement_design.py`, followed by
`close_reference_context_measurement_design.py` with their recorded pinned
plans and original launch handles. These three stages completed; do not
restart them or overwrite their outputs. Exact source counts can be exported
with `publish_reference_context_measurement_work_counts.py`; table/index paths
are exclusive outputs. Software fixtures passed 12 contexts/44 logical sides
and rejected eight rehashed false exports. Fixture native/source/journal
proofs are synthetic, not production qualification or a pilot. Resource plans
predated launch: two CPUs/16 GiB/no swap per serial stage, four GiB output
allowance and 100 GiB disk reserve; no GPU or paid resource.

This source-only bridge prepares the full phylogenetic context projection;
coverage-qualified cohorts and effects are not yet produced. Full measurement
union/coverage, common-residue and sequence-locked comparisons, domain/PAE/
prediction uncertainty, biological orthology/phylogeny and calibrated
asymmetry remain outstanding. All eight scientific aims remain incomplete.

## Full source-context coverage projection queued — September 30

`project_full_reference_context_coverage.py` joins the complete fixed source
design to the closed full-reference order/original-length coverage matrix.
The entire source record remains in `source_design`, including native guide/
context fields, every reference tie, lexical choices, parent exclusions and
duplicate model/version roles. Every reference-side/mask state retains its
source pair key, both native/usable order statuses and numerical exclusions;
unmeasured order fields are null. Missing reference models, identical models
and outside-design pairs remain separate categorical exclusions. A false
coverage flag represents eligibility, not a zero structural distance.

Five `context_policy_flags` are stored in the receipt's fixed
`context_policy_flag_order`:

1. `lexical_pair_coverage`: both lexical reference-side comparisons pass.
2. `lexical_pair_coverage_native_own`: coverage plus assignment to both
   duplicates under the source guide.
3. `lexical_pair_coverage_native_both_guides`: coverage plus assignments to
   both duplicates under both native guides.
4. `any_tie_pair_coverage_native_both_guides`: at least one original tied
   reference meets coverage and both-guide assignment.
5. `all_ties_pair_coverage_native_both_guides`: every original tied reference
   meets those requirements, with a nonempty reference set.

**Every policy also requires original parent eligibility.** Both comparisons
use the complete upstream both-order/original-length screening. Any/all-tie
policies are diagnostics; they never reselect the lexical reference. Missing
lexical structures and parent-excluded physically measured overlaps stay
explicit. These policies assess reference-side eligibility; the primary
duplicate AB comparison and independent/common-residue triads remain to be
qualified before duplication/asymmetry inference.

The full grid has 283,409 source contexts, 566,818 context/design and 1,133,636
context/design/mask records. Six screens and five policies represent 6,801,816
context/screen states and 34,009,080 decisions, plus 5,147,064 side/screen
decisions. The compressed nested full-context export remains outside Git at
`results/orthology/full-reference-context-coverage-20260930-v1`; no source
context is removed. The 240-cell `context_screen_counts.tsv` retains all four
source/parent/lexical-gene/lexical-model denominators beside passed contexts.
Counts overlap across guides/designs/masks/screens/policies.

`readback_full_reference_context_coverage.py` independently checks every full
source field, endpoint, measured/unmeasured status, numerical flag, pass and
exclusion using SQL pair lookups; it derives contextual policies through SQL
boolean/count aggregates. Both full ordinal/design/mask/screen universes, all
work counts, summary values and denominators are checked. It shares source
I/O/proof lineage only, not producer projection/policy logic. Full software
fixtures passed all states and eight rehashed false exports were rejected,
including promoted excluded parents, favorable-order passes, lexical tie
replacement, missing references, vacuous all-tie passes, cleared numerical
flags, swapped side pairs and changed denominators. Fixture closed source/
coverage/proof/journal records are synthetic, not real-data qualification.

Actual completed context-source preflight verified all 40 hash bindings and
matching current reference-inventory/primary-queue paths. Production, SQL
reader and original-journal closure are queued behind complete full-reference
coverage; no contextual production result is claimed. Plans and captured
handles are in `metadata/full_reference_context_coverage_pipeline_queued_20260930.json`.
Use its pinned wait plans and do not restart live stages or overwrite outputs.
Resources were planned before launch: two CPUs/16 GiB/no swap per serial
stage, eight GiB output allowance and 100 GiB disk reserve; 0.1–12 hours after
prerequisites is uncalibrated planning, not an ETA. No GPU or paid resource.

Full primary AB/common-triad/sequence-locked/domain/PAE/prediction controls,
biological orthology/phylogeny and calibrated inference remain outstanding.
All eight scientific aims remain incomplete; GPU predictions stay paused.

### Full three-way source-work design completed — September 30

Primary AB work and model-identity checks now accompany AR/BR for every
original context/tie. All 283,409 contexts, 214,461 ties, 428,922 logical sides
and 121,490 original availability-side links remain explicit. Independent
SQL reconstruction passed all fields, current endpoints, desired directions,
logical/physical occurrence counts and full source grids. Closure verified
59 hashes and both exact original producer/reader journals.

The full design has 31,235 ordered physical model triples, including 27,056
with source-ready logical links. Both masks and all eight AB/AR/BR order
combinations define 432,896 future correspondence states. Source readiness
requires original parent eligibility, queued AB, all three designated pair
catalogs and three distinct versioned models/model IDs; it does not establish
measured coverage or a duplication effect. Missing lexical references,
alternative ties, empty contexts and identical/outside-design models remain
unchanged. The [full workflow and evidence](full-reference-triad-work-design-20260930.md)
explain identity checks, overlapping units, native-assignment sensitivities and
remaining three-way correspondence/fitting requirements. All eight scientific
aims remain incomplete; GPU prediction stays paused.
