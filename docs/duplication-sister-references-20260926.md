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

## Independent reference-choice and path readback launched

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

The full run is active under
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
A successful final receipt is required; launch alone is not a passed audit.
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
