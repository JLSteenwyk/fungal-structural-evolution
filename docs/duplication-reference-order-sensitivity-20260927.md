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
