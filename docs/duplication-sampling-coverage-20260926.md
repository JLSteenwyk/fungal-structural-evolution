# Structural coverage of terminal duplication candidates

`scripts/summarize_duplication_sampling_coverage.py` independently reconstructs
all reported Terminal events with exactly one gene on each side from both
complete native duplication tables. It joins both genes to the same frozen
AlphaFold structural bridge used by the structural-comparison queues and
compares event identity, singleton status and model-presence counts against
the original coverage exports. This includes events with no models or only one
model, which do not enter the current structural pair queue.

The output contains an event-disposition table and summaries by taxon and
family. Every candidate-manifest taxon is retained, including zero qualifying events.
Zero denominators produce blank fractions, not zero structural coverage.
Counts distinguish neither/one/both genes modeled and identical-model cases.
Species names, study roles and lineage strings come from the pinned sampling
manifest. Profiles and MAFFT are kept separate, not treated as independent
replicates. Reported terminal singleton-side events are a restricted candidate
class, not the full biological duplication history.

The immediate purpose is to quantify ascertainment of the modeled candidate
set before downstream divergence/asymmetry interpretation. Missing models are
not biological absences. Counts do not establish missing-at-random sampling,
correct ascertainment bias or account for phylogenetic dependence. The frozen
bridge excludes later retrievals and ESMFold structures; the table does not
claim current total structure availability. It also does not independently
validate every missing-model candidate's gene-tree node.

Plan: `metadata/duplication_sampling_coverage_plan_20260926.json`.
Launch: `metadata/duplication_sampling_coverage_launch_20260926.json`.
Output: `results/orthology/duplication-sampling-coverage-20260926-v1/`.
Run with:

```bash
python scripts/summarize_duplication_sampling_coverage.py \
  --plan metadata/duplication_sampling_coverage_plan_20260926.json
```

The job completed. Resources were one CPU, 8 GiB RAM, no swap and an estimated
1 GiB output. The 0.1–4 hour planning range is uncalibrated. No GPU or paid
resources are used. Independent aggregate readback passed. Verified results:

| Guide | Terminal singleton-side events | Neither model | One model | Both models | Taxa with both models |
|---|---:|---:|---:|---:|---:|
| Profile | 467,663 | 345,890 | 12,528 | 109,245 | 153 |
| MAFFT | 467,690 | 345,933 | 12,529 | 109,228 | 153 |

All 526 reconciled taxa have reported events in this restricted class, but only 153
contribute two-model events. Two-model coverage is approximately 23.4% in
each guide. Both include 6,354 same-model events. Profile and MAFFT respectively
have 48,661 and 48,711 families with events; only 20,492 and 20,474 have
two-model events. This substantial ascertainment limits generalization across
lineages and families. Counts are archived in
`metadata/duplication_sampling_coverage_completed_20260926.json`; all source
and output hashes were rechecked after completion.


## Aggregate readback and manifest-versus-analysis membership

The independent grouped-event readback passed all 935,353 event records,
1,054 taxon/guide rows and 97,372 family/guide rows. Every model-presence class,
same-model flag, group count, fraction and taxon metadata field matched.
It uses exported event data rather than the producer's streaming counters;
it does not repeat the producer's native-event/bridge joins.

The raw candidate manifest contains **527 entries**, while both audited
reconciliation species lists contain **526 taxa**. The extra entry is
*Saccharomyces jurei* (F1987369), explicitly excluded because no usable annotated
proteome was acquired. It remains in the candidate manifest for possible
recovery. Thus the original coverage table has one blank-denominator row per
guide that means outside the analysis, not a reconciled species with biological
absence of duplications. The producer receipt's `sampled_taxa=527` refers to
candidate-manifest rows, not the analyzed reconciliation cohort.

A separate immutable table now adds `in_reconciliation`, `denominator_status`
and `exclusion_reason`, preserving every original field:
`results/orthology/duplication-coverage-membership-20260926-v1/taxon_coverage_with_membership.tsv`.
Both guides have 526 `events_observed` rows and one `not_in_reconciliation`
row. Membership was checked against each audited native SpeciesIDs file and
`metadata/sampling_exclusions.json`. Use this explicit membership table for
future plots and downstream ascertainment summaries.

Reproduce with:

```bash
python scripts/readback_duplication_sampling_coverage.py --source-plan metadata/duplication_sampling_coverage_plan_20260926.json --output <fresh-readback.json>
python scripts/annotate_duplication_coverage_membership.py --coverage results/orthology/duplication-sampling-coverage-20260926-v1 --aggregate-readback <fresh-readback.json> --exclusions metadata/sampling_exclusions.json --profile-tree-audit metadata/expanded_resolved_tree_profile_completed_readback_20260923.json --mafft-tree-audit metadata/expanded_resolved_tree_mafft_completed_readback_20260923.json --output <fresh-membership-directory>
```

Completed readback and membership receipts are archived in
`metadata/duplication_sampling_coverage_readback_completed_20260926.json`.


## Lineage coverage figure

[Coverage figure](figures/duplication_lineage_coverage_20260926.png)
([PDF](figures/duplication_lineage_coverage_20260926.pdf),
[SVG](figures/duplication_lineage_coverage_20260926.svg)) summarizes 27 broad
manifest lineage groups separately for both guides. The left panel uses events
as the denominator; the right uses reconciled taxa. The excluded candidate
*S. jurei* is omitted. All 54 rows and 108 percentages were independently
reconstructed from the membership table and output hashes checked; the final
render was visually inspected. Verification is recorded in
`metadata/duplication_lineage_coverage_figure_completed_20260926.json`.

The different denominators matter: Blastocladiomycota has approximately 57%
event coverage but only two of nine taxa contribute modeled pairs. Coverage
within large lineages is also uneven: only 56 of 234 Ascomycota contribute
pairs. Four ingroup lineages have no two-model events in this frozen candidate
set. These observations describe structural ascertainment of a restricted
event class, not biological absence or current total model availability.

Reproduce in fresh output locations:

```bash
python scripts/plot_duplication_sampling_coverage.py \
  --membership results/orthology/duplication-coverage-membership-20260926-v1 \
  --output <fresh-output-directory> \
  --figure-prefix <fresh-figure-prefix>
```

The figure script produces a numeric `lineage_coverage.tsv`, provenance receipt,
and PNG/PDF/SVG figures. No new prediction or phylogenetic correction is performed.

## September 28 expanded-catalog join

The fixed 935,353 terminal singleton-side event records have now been joined
to the September 28 catalog of 1,955,694 protein/model links. This separate
output preserves the original event identities and old model assignments:
`results/orthology/duplication-expanded-coverage-20260928-v1/`.
The producer and independent SQL readback both finished successfully. Every
original event field, new model assignment, coverage class and aggregate count
was checked across all 935,353 rows.

Verified counts are 141,724 two-model events for the profile guide
and 141,685 for MAFFT, versus 109,245 and 109,228 previously. Both guides now
have two-model candidates in 210 taxa, versus 153 previously. No previously
two-model event lost that coverage. These are source-availability counts, not
confidence-qualified comparisons or evidence of structural divergence.
The old figures and comparison queues retain their original frozen scope.

Reproduce in a fresh output directory using
`scripts/refresh_duplication_model_coverage.py --plan
metadata/duplication_expanded_coverage_plan_20260928.json`.
The plan pins source tables, receipts and script; the launch record is
`metadata/duplication_expanded_coverage_launch_20260928.json`. The bounded job
used one CPU, 4 GiB RAM, no swap, no GPU and no paid resources, with a 1 GiB
output allowance and an uncalibrated 0.02–1 hour planning range.

The independent verifier is `scripts/readback_expanded_duplication_coverage.py`;
run it with the same `--plan` and a fresh `--output` JSON path. It reconstructs
both protein/model joins in SQLite instead of using the producer's lookup and
recomputes every summary with SQL. Completion evidence and source hashes are
in `metadata/duplication_expanded_coverage_completed_20260928.json`.
Two-model coverage is 30.30% for profile and 30.29% for MAFFT. Both contain
7,294 events whose two proteins map to the same model; model availability
does not imply two independent structural predictions. Readback relies on
the earlier audited native-event selection and sequence catalog and does not
repeat coordinate or biological-event validation.

All 1,054 taxon/guide membership rows have also been updated in
`results/orthology/duplication-expanded-membership-20260928-v1/`.
Every original event denominator and cohort label is preserved; the excluded
*S. jurei* row retains its blank coverage fraction. Streaming event counts were
checked against SQL aggregation, and every serialized output row was read back.
Both guides gain 57 taxa with at least one two-model event: 32 Ascomycota,
12 Basidiomycota, five Mortierellomycota, five Glomeromycota, one Olpidiomycota,
one Mucoromycota and one Amoebozoa outgroup. These counts describe newly
available candidates, not independent ecological transitions or qualified
structural-divergence evidence.

Reproduce with `scripts/summarize_expanded_duplication_taxa.py --readback
results/orthology/duplication-expanded-coverage-readback-20260928-v1.json
--membership results/orthology/duplication-coverage-membership-20260926-v1
--output <fresh-directory>`. The bounded job used one CPU, 2 GiB RAM and no
swap or GPU, finishing successfully in about 10 CPU seconds. Source/output
checksums are recorded in `metadata/duplication_expanded_taxa_completed_20260928.json`.

The updated [lineage figure](figures/duplication_lineage_coverage_20260928.png)
([PDF](figures/duplication_lineage_coverage_20260928.pdf),
[SVG](figures/duplication_lineage_coverage_20260928.svg)) uses this expanded
membership table. All 54 guide/lineage rows and 108 percentages were checked
with a separate counter-based aggregation, and the PNG was visually inspected.
Evidence is in `metadata/duplication_expanded_lineage_figure_completed_20260928.json`.
Reproduce with the existing `scripts/plot_duplication_sampling_coverage.py`,
passing the expanded membership directory and fresh output/figure paths.

The two denominators remain different: about 42.1% of Ascomycota candidate
events have both models, but only 88 of its 234 taxa contribute covered pairs.
Olpidiomycota now has nearly complete event coverage, but represents one taxon.
Aphelidiomycota, Sanchytriomycota and Calcarisporiellomycota still have no
two-model candidates in this frozen event class; that is not evidence of
biological absence or a survey of every available structure.

### Expanded comparison resource inventory

The September 28 inventory contains 134,812 distinct model/version pairs
across 268,821 event links, plus 14,588 same-model event links retained in
the source table. Of these distinct pairs, 103,199 also occur in the old
103,200-pair queue and 31,613 are new. Comparing all new pairs in both input
orders would require 63,226 alignments per mask, or 126,452 with two masks,
before eligibility filtering. These are counts, not a measured runtime ETA.

The one old pair absent from the expanded set involves F13290 proteins
CAD7067269.1 and CAD7067533.1 in both guides. The latter's selected model
changed from AF-A0A177T084-F1 to AF-A0A9N8M122-F1; two-model event coverage
was retained. Matching model IDs and versions alone does not justify reusing
old alignment results: coordinate bytes, masking and settings must also match.

Reproduce with `scripts/estimate_expanded_duplication_pairs.py --readback
results/orthology/duplication-expanded-coverage-readback-20260928-v1.json
--old-queue-evidence metadata/duplication_model_pair_queue_completed_20260926.json
--output <fresh-directory>`. All 935,353 event rows were scanned; every
distinct-pair multiplicity was checked against SQL aggregation and every
serialized pair row was checked. Source hashes were verified before and
after execution. The one-CPU, 4-GiB, zero-swap job finished successfully in
18 CPU seconds, with a 566-MB peak. Its planning allowance was 100 MiB output
and 0.01–1 hour; no GPU, paid resource or new alignment job was used.
Evidence is in `metadata/duplication_expanded_pair_estimate_completed_20260928.json`;
the pair table and receipt are in
`results/orthology/duplication-expanded-pair-estimate-20260928-v1/`.

This inventory is not an executable comparison queue. Expanded candidate
tree-node correspondence, coordinate/confidence checks and masking eligibility
remain prerequisites. Source event associations remain intact, and deduplicating
model pairs for computation does not create independent biological observations.

### Expanded native-event inventory and tree review

The unchanged `inventory_duplication_structure_coverage.py` has completed a
new full native-event inventory using the audited September 28 family bridge.
The output at `results/orthology/duplication-structure-coverage-20260928-v1/`
retains all 1,204,638 profile and 1,204,919 MAFFT events, including complex and
uncovered events. Its terminal two-model exports contain 141,724 and 141,685
candidates, respectively, each spanning 210 taxa. Source-plan binding, artifact
hashes and successful service termination are archived in
`metadata/duplication_structure_coverage_completed_20260928.json`.

Every candidate identity, taxon, gene/model assignment and same-model flag was
compared against the independent expanded catalog join, scanning all 935,353
terminal singleton-side records. The initial comparison exposed a formatting
difference: native gene order versus sorted gene order. The corrected checker
sorts paired gene/model assignments together, preserving the mapping. Fixtures
verify order invariance and rejection of changed assignments, genes, taxa and
flags. All 283,409 candidates passed. This does not check every noncandidate
coverage field or establish structural eligibility.

Reproduce the inventory with `scripts/inventory_duplication_structure_coverage.py
--plan metadata/duplication_structure_coverage_plan_20260928.json` in a fresh
output location. The plan allows one CPU, 8 GiB RAM, zero swap, 2 GiB output and
0.1–4 hours; the run completed in 86 CPU seconds with an 804-MB peak. Reproduce
the candidate comparison with `scripts/check_expanded_duplication_candidate_export.py
--coverage results/orthology/duplication-structure-coverage-20260928-v1
--expanded-readback results/orthology/duplication-expanded-coverage-readback-20260928-v1.json
--output <fresh-json-path>`; run its focused identity fixtures with
`scripts/check_duplication_candidate_export_identity.py`.

The full tree review completed successfully as
`fungal-duplication-candidate-tree-review-20260928.service`, using unchanged
`validate_duplication_structure_candidates.py` with
`metadata/duplication_candidate_tree_review_plan_20260928.json`.
It independently rebuilds eligibility from native events and the frozen model
bridge, then checks reported nodes in 50,063 guide-specific candidate-bearing
family trees. Both descendant identity and direct-tip-child status are retained;
missing or mismatching nodes remain explicit. The plan caps one CPU and 8 GiB
RAM with zero swap, allows 1 GiB output and an uncalibrated 0.1–6-hour runtime.
Launch identity is in `metadata/duplication_candidate_tree_review_launch_20260928.json`.
All 141,724 profile and 141,685 MAFFT candidates match their exact reported
nodes and have the expected two direct tip children. Both guides share 141,302
protein pairs, with 422 profile-only and 383 MAFFT-only pairs. Every exported
source candidate field was checked, and all summary counts and guide-pair
membership rows were recomputed after successful termination. Evidence is in
`metadata/duplication_candidate_tree_review_completed_20260928.json`.
No GPU work or structural alignment was launched; rooting, support, duplication
biology and structural effects remain separate requirements.

The unchanged `prepare_duplication_structure_pairs.py` completed the
expanded model/version queue at
`results/structural_comparisons/duplication-model-pair-queue-20260928-v1/`.
Its complete command, source pins and resource allowance are recorded in
`metadata/duplication_model_pair_queue_plan_20260928.json`; launch identity is
in `metadata/duplication_model_pair_queue_launch_20260928.json`.
This CPU-only job retains every reviewed event association and explicitly
flags same-model rows. It is limited to one CPU, 8 GiB RAM, zero swap, with
1 GiB output and an uncalibrated 0.02–2-hour planning allowance. It finished in
54 CPU seconds. The queue contains 283,409 event links, 14,588 same-model links,
276,682 total models and 134,812 distinct model/version pairs. The 269,393
models used in distinct pairs occupy 94,489,735,215 coordinate bytes by file
size. This is a storage measurement, not raw-content validation.

The separate `readback_duplication_model_pair_queue.py` completed successfully
in 47 CPU seconds with a 1.4-GiB peak. It checked every original reviewed
event field, pair key/status and membership, every complete frozen catalog
record and all corresponding counts. Sources were hashed before and after;
all output rows were checked. It does not independently repeat the original
protein-to-model join or raw coordinate parsing. Reproduce with `--queue`,
`--review`, `--catalog` pointing to the above September 28 outputs and a fresh
`--output` JSON path. Completion evidence is in
`metadata/duplication_model_pair_queue_completed_20260928.json`.

Raw-coordinate validation is now running for all 276,682 models, including
same-model cases, through unchanged `validate_duplication_coordinates.py` and
`metadata/duplication_coordinate_validation_plan_20260928.json`. The source
files total 96,596,573,260 bytes; maximum sequence length is 1,567 residues.
The job checks raw byte hashes, sequence/atom identity, C-alpha coordinates and
per-residue confidence, retaining content rejections explicitly. It uses four
CPU workers, 16 GiB RAM and zero swap, with 1,000-model checkpoints, 32 GiB
output allowance and a 100-GiB free-disk reserve. The 0.5–24-hour planning range
is uncalibrated. Launch identity is in
`metadata/duplication_coordinate_validation_launch_20260928.json`; output is
`results/structural_comparisons/duplication-coordinate-validation-20260928-v1/`.
Coordinate validation, independent numeric readback and structural comparisons
remain pending. No GPUs or paid resources are used. Matching model IDs/versions
alone still does not justify reuse of previous coordinate-derived results.
