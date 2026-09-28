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
