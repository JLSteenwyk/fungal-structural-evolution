# Reference comparison sensitivity to input order

Both alignment directions are now compared across every reference pair/mask
with two numerically usable mappings: 32,541 full-mask and 30,346 pLDDT70 pairs.
All 125,774 native residue mappings were reconstructed from hashed checkpoints,
with reverse inputs converted back to the same endpoint orientation. Four
pLDDT70 pairs with only one usable direction and 2,191 with neither remain in
the parent summary; they cannot contribute a two-direction contrast.

| Cohort | Pairs | Same residue correspondence | Different correspondence |
| --- | ---: | ---: | ---: |
| Full, all usable pairs | 32,541 | 32,518 | 23 |
| Full, shared mask cohort | 30,346 | 30,331 | 15 |
| pLDDT70, shared mask cohort | 30,346 | 30,287 | 59 |

The shared cohort uses exactly the same pairs in both masks. Of the full-mask
23 differences, 20 change the number of aligned pairs and three retain the same
count but change correspondence. For pLDDT70, these counts are 27 and 32.
Equal aligned length alone therefore does not establish equivalent residue
mapping. All identical correspondences yield recomputed RMSDs agreeing within
1e-10 Å; most observed order differences are numerical rounding.

The largest absolute RMSD difference between directions is 1.0165 Å for full
inputs and 2.6080 Å for pLDDT70 inputs. The largest endpoint-A-normalized native
TM-score differences are 0.03709 and 0.52318 respectively. These extremes concern
potentially different aligned residues. They are not structural uncertainty
bounds or biological changes, and the mask comparison does not isolate a causal
confidence effect. No direction was selected as preferred, and no average
replaced the two values.

## Verification and reproduction

Run `scripts/summarize_reference_order_sensitivity.py`, followed by
`scripts/readback_reference_order_sensitivity.py`. Existing output directories
and proof files are immutable; use new paths for a repeat execution. Full output
is `results/structural_comparisons/duplication-reference-order-sensitivity-20260927-v1`.
The checker independently rebuilt all 125,774 mappings using cumulative-index
arrays and checked all 62,887 pair/mask rows, shared-cohort membership, metric
differences and 27 quantile rows (sorted interpolation versus NumPy quantiles).

Versioned evidence is in
`metadata/duplication_reference_order_sensitivity_{receipt,readback}_20260927.json`
and `metadata/duplication_reference_order_sensitivity_quantiles_20260927.tsv`.
The full cohort still needs scientific coverage/confidence qualification and
integration with primary duplicate/background comparisons. Input-order stability
alone does not establish prediction accuracy or evolutionary asymmetry.

## Full expanded reference order and coverage workflow queued — September 30

The earlier results above concern the September 27 cohort. The full expanded
reference design now contains 55,701 pairs and 222,804 mask/order states. New
scripts `screen_full_reference_orders_and_coverage.py` and
`readback_full_reference_orders_and_coverage.py` consume the source/journal-
closed full measurement union. They normalize usable metrics to current model
endpoints while retaining original source orders, checkpoint hashes, all
exclusions, raw aligned lengths and exact union-row hashes. Both usable orders
are required for each of the six existing original-protein coverage screens;
confidence masking does not change the denominator. Retained-input coverage
is separately labeled. Neither a preferred order nor an average is used.

Software fixtures passed 20 states and independently reconstructed all fields
and 60 decisions, including reversed source directions, unequal versions,
nonlexical current endpoints, denominator and decimal-boundary cases, all
native statuses and all numerical exclusions. Seven rehashed false exports
and one incorrect original-length source were rejected. Fixture source/proof/
journal records are explicitly synthetic and do not qualify real results.

Actual complete-source preflight passed all 83,207 models and both frozen input
manifests (602,972 states). Evidence is
`metadata/full_reference_order_coverage_input_preflight_verified_20260930.json`:
60 source/artifact bindings and the original process completion journal. Full
production, independent reconstruction and two-journal closure are queued;
see `metadata/full_reference_order_coverage_pipeline_queued_20260930.json`.
Run each stage through its pinned wait plan, using the captured original
dependency handles; do not rerun an already live stage or overwrite its output.
Resources are two CPUs/16 GiB/no swap per serial stage, no GPU or charges, an
8 GiB output allowance and 100 GiB disk reserve. The 0.1–12 hour planning
range after prerequisites is uncalibrated, not an ETA.

Production outputs will be outside Git at
`results/structural_comparisons/full-reference-order-coverage-20260930-v1`;
completion will require full independent field/threshold/count reconstruction
and both exact original process journals. This is metric/coverage screening;
it does not yet reconstruct the full native residue mappings for comparison
of correspondence between orders. The full 121,490 reference-side links and
283,409 target contexts/designs remain to be projected without dropping missing
references or parent exclusions. Common triads, sequence-locked/domain/PAE/
prediction controls, biological orthology/phylogeny and calibrated asymmetry
remain required. All scientific aims remain incomplete.
