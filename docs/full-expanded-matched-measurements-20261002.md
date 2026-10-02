# Complete expanded matched measurement handoff

The expanded dataset now has a full case index and a validated measurement
join. This stage connects the measured protein comparisons to the original
gene contexts needed for sequence–structure and duplication models. It does
not yet estimate a calibrated duplication effect or correct for phylogenetic
dependence.

The case producer and independent SQL reader passed all **4,250,692 fixed
selections**, **1,133,636 target/policy records**, and **56,965,652 unmatched
decisions**. Final [provenance closure](../metadata/full_matching_case_index_completed_20261002.json)
passed 1,830,493 source/artifact bindings and both original completion/resource
journals. Their complete census contains **75,188 logical cases** and
**37,600 distinct ordered physical target/control pair combinations**.
The two guides contribute 37,588 and 37,600 logical cases; these overlapping
records are dependent. Matching scenarios and policies reuse some cases up
to 216 times. Neither the selections nor the distinct model combinations are
independent evolutionary samples.

The selected records span **3,331 families and 199 focal taxa per guide**.
They are an ascertained subset of the 526-entry project design. At the logical
case level, 70,221 have distinct models on both sides, 4,788 use the same model
on both sides, 139 use the same model only for the background, and 40 only for
the target. The latter cases remain explicit and receive no invented zero
distance. All unmatched scenarios remain in their original closed status
table; quality filtering never causes rematching.

## Sources, identity and outputs

The [data dictionary](../metadata/expanded_matched_measurement_data_dictionary_20261002.tsv)
lists every concrete field in the case index, selection membership, directed catalog and contrast
table, including types, units, foreign keys, quality bits and missing-value
semantics. It is reproducible with `write_expanded_measurement_dictionary.py`.

| Artifact | Full scope and interpretation |
| --- | --- |
| `case_index.tsv.gz` | 75,188 target/background node combinations, retaining guide, families, focal taxon, gene-tree node, all four genes, background taxa, physical pair identities, sequence distances, reuse and quality bits |
| `selection_case_links.tsv.gz` | All 4,250,692 original selections, including exact score, ties, scenario, policy and endpoint mapping, with source ordinal and logical/physical case references |
| Original target/policy coverage status | All 1,133,636 records and 56,965,652 unmatched decisions; referenced with its original checksum, fully replayed for membership |
| `directed_measurements.tsv.gz` | Complete 1,123,936 states: 539,248 target and 584,688 background states, each physical pair/mask/order, with original checkpoints and sources |
| `case_mask_contrasts.tsv.gz` | Complete 150,376 case/mask rows and 1,203,008 potential order-pair/outcome cells; full and pLDDT70 are physical masks, while their intersection is an eligibility gate |

`case_id` hashes the ordered target/background **node** identifiers with the
`fixed-matched-logical-case-v1` namespace. `physical_case_id` hashes the ordered
target/background **model-pair** identifiers with a separate namespace.
Genes sharing models stay distinct. Case identity never depends on favorable
quality or observed structural outcomes. Each future case/mask `row_identity`
hashes its case ID and mask with `expanded-matched-measurement-row-v1` for
covariance and model-input linkage. Both source graph files remain part of the
full closed provenance, including endpoint model versions, confidence, lengths,
sequence hashes and coordinate hashes.

The directed catalog preserves original native order labels separately from
the canonical input direction. It retains raw numeric diagnostics for excluded
native-aligned states with **`numerical_usable=0`**. Native failures keep blank
metrics; actual usable zero RMSD stays zero. Target coverage intentionally
blanks excluded measurements, whereas background coverage retains raw aligned
lengths with explicit exclusion flags. The catalog respects both original
contracts. Its raw excluded numbers must never be used as valid model outcomes.

`coverage_left/right` are the native alignment-input fractions. The separate
`original_coverage_a/b` use the complete protein lengths. `order=0` means the
canonical model A is native left and B is right; order 1 reverses those inputs.
`source_order` is the original checkpoint's label and can differ from `order`.
Native TM metrics remain labelled as native scores, not independently optimized
TM scores. Checkpoint keys map every normalized state back to the exact original
result and versioned input chain.

## Contrasts and numerical validation

Every logical case retains both target orders and both control orders. For
each mask, the join calculates all four target-minus-control comparisons for
independently checked RMSD and native TM dissimilarity, defined as
`1 - (TM_left_native + TM_right_native)/2`. The endpoint average is symmetric;
its protein-length normalization and prediction-source limitations remain.

Means require both usable orders on the corresponding side. The mean
target-minus-background contrast requires both complete means. A complete
minimum/maximum/span envelope requires all four target/control order cells.
Partial available deltas remain visible, but cannot become a complete envelope;
failed and same-model states remain null. These envelopes describe dependent
alignment-order sensitivity, not confidence intervals, ancestral polarity,
significance or calibrated evolutionary effects.

The case index reader independently reconstructs every membership and case
aggregate in SQLite, including scenario 54 and unmatched partitions. The
catalog reader independently reconstructs every source, order, numeric value,
numerical exclusion and quality projection. The join reader separately
implements state gating, unit conversion, averaging and all four differences
using 50-digit Decimal arithmetic, allowing at most `2e-14` times the relevant
input scale for ordinary float rounding. Full source/artifact hashes and two
original invocation-linked completion/resource journals gate each final
handoff. Producer receipts or collected systemd success defaults alone do not
establish completion.

Software contracts passed 23 invalid case-index exports, 20 invalid catalog
exports, and 19 invalid join exports, with full replay after interruption and
refusal to restart a completed stage. Fixtures test the new software; prior
native fits, matching and journal proofs in fixtures are synthetic contracts,
not a biological pilot or production acceptance.

## Execution and reproducibility

The [case index plan](../metadata/full_matching_case_index_plan_20261002.json)
has a successful full producer/reader and completed original provenance closure. The
[corrected catalog v3 plan](../metadata/full_expanded_measurement_catalog_plan_20261002_v3.json)
and [join v2 plan](../metadata/full_expanded_case_measurements_plan_20261002_v2.json)
have completed directed-catalog provenance closure for all 1,123,936 states,
with 2,369,743 bound hashes and both original journals. The full join passed
[separate Decimal readback and final closure](../metadata/full_expanded_case_measurements_v2_completed_20261002.json)
for all 150,376 case/mask rows, with 2,369,783 bound hashes and both original
journals. Complete four-order envelopes exist for 70,221 full-mask and
70,015 pLDDT70 cases per outcome; 9,934 same-model case/mask rows retain null
measurements. Exact identities remain recorded in the versioned launch
inventories. Existing native fits and other scientific jobs were not changed.

Resources were recorded before each launch: two CPU equivalents, 32 GiB memory,
no swap, one BLAS thread, 100 GiB free-space reserve and 16–32 GiB output/scratch
allowances. The 1–12 hour planning allowance per stage is uncalibrated and
excludes dependency waiting; it is not a project completion ETA. There are no
new GPU jobs or paid resources. Large measurements and complete proof maps
remain under `results/`, outside Git.

The first catalog attempt stopped because its downstream matching archive
did not list the original target manifest. The corrected source loader expands
the original verified usable-summary receipt chain through native/diagnostic/
geometry receipts, the complete manifest and all 539,248 native checkpoints.
The second attempt exposed the target coverage table's intentional blanking
of excluded values; v3 preserves that blanking while retaining quarantined raw
diagnostics. Both failed attempts and the original dependent join failures
remain recorded; their plans/scripts/results/journals were not overwritten.
See the [first attempt](../metadata/full_expanded_measurement_catalog_attempt_20261002.json)
and [second attempt](../metadata/full_expanded_measurement_catalog_attempt_20261002_v2.json).

Reproduce each stage with its producer and reader and `--plan` pointing to the
versioned plan, using a new plan/output identity for a new run. The original
plan pins its scripts, fixtures, resource estimate and source handoffs. Do not
edit a launched plan or rerun its completed output. Final completion locators
are written by the original queued closure jobs after the full reader and
journal/hash checks pass.

The [complete expanded covariance index/factor workflow](full-expanded-covariance-20261002.md)
has completed full source, entity, rank and raw-tree covariance checks for this
entire cohort. Next required work includes the expanded model-input design; full mask/order/screen/matching/reuse sensitivities;
accepted species/gene/reconciliation uncertainty; dependence and ascertainment
modelling; predictor, domain, PAE and orientation controls; calibration and
multiple-testing treatment. Older 375,350 whole-protein fits retain their
earlier frozen cohort and are not proof of fitting this expanded dataset.
All eight biological aims remain incomplete.
