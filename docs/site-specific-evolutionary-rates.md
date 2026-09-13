# Conditional site-specific evolutionary rates

Completed 288 fixed-topology fits across all 72 local-source paired markers,
exporting 66,284 relative site-rate estimates. The models are AA LG+F+G4,
3Di AF+G4, AF+F+G4 and LLM+G4. Each uses the original paired observations and
fixed AA topology; branch lengths and model parameters are reestimated.

IQ-TREE's [empirical-Bayes rate export](https://www.iqtree.org/doc/Advanced-Tutorial#inferring-site-specific-rates)
is the posterior mean of discrete rate categories. The assigned modal category
and its rate are also retained. These multipliers describe relative rates
within each fitted model, not substitutions per year, physical displacement,
positive selection or calibrated uncertainty. Four-category Gamma is a baseline
specification here; model adequacy and alternative heterogeneity models have
not yet been established. Comparisons must account for shrinkage and different
alphabets rather than equating raw AA and 3Di rate values.

## Execution and validation

The first attempt stopped at the parser: installed IQ-TREE 3 writes `Cat` and
`C_Rate`, while the documentation example uses `Category` and `Categorized_rate`.
The parser now supports both schemas; tests also reject malformed site grids,
nonfinite/negative rates and out-of-range categories. The failed v1 directory
is preserved without a completion receipt; all results below use v2.

The full audit checked all fit/source hashes, 288 model/topology identities,
66,284 rate/category entries, Gamma-category bounds and site-likelihood sums.
Every likelihood sum agrees with its report within the precision implied by
its printed numbers. The maximum absolute sum discrepancy is 0.00567 log units.
Compared with previous point fits, refitted log likelihoods differ by -0.0017
to +0.0363. These small differences are recorded rather than treated as exact
parameter replay. All 248 fits with warnings remain explicit, including
near-zero branch and missingness/state-frequency warnings.

The audit does not independently implement the likelihood calculation or
recompute posterior rate weights. Prediction artifacts, topology/model
uncertainty, source effects and site dependence remain relevant. Joining rates
to the exposure/functional-site annotations and performing controlled comparative
analyses is still pending.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/estimate_paired_site_rates.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --audit results/phylogeny/paired-fit-audit-esmfold-v1 --output results/phylogeny/paired-site-rates-esmfold-v2
OPENBLAS_NUM_THREADS=1 python scripts/audit_paired_site_rates.py --rates results/phylogeny/paired-site-rates-esmfold-v2 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --output results/phylogeny/paired-site-rates-audit-esmfold-v1
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_site_rate_parser.py
```

Full site tables and individual fit artifacts are outside Git. Provenance,
configuration, fit summaries and warnings are in `metadata/esmfold_site_rate_*`.
Resource planning used four single-thread fits, 8 GB RAM and 5 GB disk headroom
on the existing authorized host; no new paid resources were used.


## Matched analysis table completed

`scripts/assemble_site_evolution_frame.py` now combines all 16,571 sites and
66,284 rate estimates with their existing exposure/parsimony/topology-sensitivity
fields. It adds observed taxon fractions, focal CA confidence minima/medians/
quartiles, each fitted Gamma shape and fit warning count. All 714,936 original
taxon-site observations are accounted for. Site identities and copy-review flags
are retained; no filtering, regression or significance testing is performed.

The assembly requires the completed rate and topology audits and checks common
input/fit lineage. Full readback verified every inherited field and rate/category
value. Confidence statistics were recalculated with sorted interpolation and
standard-library medians. Every model parameter/warning covariate, tree taxon
count and marker summary matched its source. Thus the next statistical step can
use matched measurements without silently changing the site population.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/assemble_site_evolution_frame.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --topologies results/phylogeny/site-parsimony-topologies-esmfold-v1 --topology-audit results/phylogeny/site-parsimony-topologies-audit-esmfold-v1 --diagnostic results/phylogeny/site-parsimony-exposure-esmfold-v1 --rates results/phylogeny/paired-site-rates-esmfold-v2 --rate-audit results/phylogeny/paired-site-rates-audit-esmfold-v1 --projection results/structural_annotations/paired-accessibility-esmfold-v1 --output results/phylogeny/site-evolution-frame-esmfold-v1
```

The full table is outside Git. Its checksum, provenance, readback and coverage
summaries are in `metadata/esmfold_site_evolution_frame_*`. Availability of these
covariates does not mean their effects have been adjusted for. Estimated rates
and exposure remain dependent and conditional on prediction/model choices.

See [paired-marker model diagnostics](paired-model-adequacy.md) for the completed unfiltered symmetry screen and its substantial limitations in test availability.


## Four-category FreeRate sensitivity: execution complete

The full 72-marker × four-model FreeRate run completed at
`results/phylogeny/paired-site-rates-freerate-esmfold-v1`. It changes only the
heterogeneity specification from `+G4` to `+R4` in the original matrix/frequency
models; all original paired observations, fixed AA topologies and seeds are
retained. Branch lengths, category rates and weights are reestimated. Resource
planning is four single-thread workers, 8 GB RAM, 5 GB disk and 0.1–4 hours on
the existing host. No paid resources are used.

[IQ-TREE's rate-heterogeneity documentation](https://iqtree.github.io/doc/Substitution-Models#rate-heterogeneity-across-sites)
defines FreeRate as relaxing the Gamma assumption. This comparison assesses
sensitivity to that assumption; it does not establish model adequacy. Increased
likelihood alone is insufficient evidence because FreeRate adds flexibility.
No likelihood-ratio p-values or automatic model selection are planned here,
given unresolved feature dependence and sparse-state diagnostics.

The rate-export script accepts `--heterogeneity R4`; the default remains G4.
The revised full auditor was rerun across all 288 completed Gamma fits. Every
previous fit-summary value matched, and the site-rate and warnings tables were
byte identical. The summary adds an explicit `rate_heterogeneity` column and
leaves Gamma alpha blank for FreeRate. Two rate-parser tests passed.

`scripts/compare_site_rate_heterogeneity.py` is prepared to compare both audited
runs on identical inputs. It exports matched site and branch values and
within-fit rank correlations, tree-length totals and likelihood changes. It
retains lower FreeRate likelihoods as optimization diagnostics. Tied ranks
were checked against SciPy; constant-vector cases remain unavailable. A negative
test correctly rejected Gamma output supplied as the FreeRate input before
creating output. The full comparison and readback are now complete; seven optimization
diagnostics remain unresolved, as described below. Category numbers are not matched across
models, and relative site multipliers are not absolute evolutionary rates.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/estimate_paired_site_rates.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --audit results/phylogeny/paired-fit-audit-esmfold-v1 --heterogeneity R4 --output results/phylogeny/paired-site-rates-freerate-esmfold-v1
# After the complete FreeRate receipt exists:
OPENBLAS_NUM_THREADS=1 python scripts/audit_paired_site_rates.py --rates results/phylogeny/paired-site-rates-freerate-esmfold-v1 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --output results/phylogeny/paired-site-rates-freerate-audit-esmfold-v1
OPENBLAS_NUM_THREADS=1 python scripts/compare_site_rate_heterogeneity.py --gamma results/phylogeny/paired-site-rates-esmfold-v2 --gamma-audit results/phylogeny/paired-site-rates-audit-esmfold-v2 --free results/phylogeny/paired-site-rates-freerate-esmfold-v1 --free-audit results/phylogeny/paired-site-rates-freerate-audit-esmfold-v1 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --output results/phylogeny/paired-rate-heterogeneity-esmfold-v1
```


### Completed comparison and optimization diagnostics

All 288 FreeRate fits completed and passed the full rate-output audit, including
66,284 rate entries, model/topology identities, category bounds and site-likelihood
accounting. Warnings are retained for 249 fits. The matched comparison contains
66,284 site pairs and 36,488 branch pairs. Full readback matched every site value
to its raw rate export and every branch value to the source trees through the
shared topology helper. All 576 within-fit rank correlations were separately
checked with SciPy; median differences, tree totals, likelihoods and flags also
matched. This does not independently implement likelihoods, posterior weights
or tree parsing.

| Matrix/frequency specification | Median site-rate rank correlation | Median branch-length rank correlation | FreeRate likelihood lower by >0.1 |
| --- | ---: | ---: | ---: |
| AA LG+F | 0.997972 | 0.999486 | 1/72 |
| 3Di AF | 0.996468 | 0.998500 | 1/72 |
| 3Di AF+F | 0.995151 | 0.997607 | 3/72 |
| 3Di LLM | 0.995608 | 0.998677 | 2/72 |

These medians describe all 72 markers equally; high rank agreement does not
establish equal absolute lengths, stable acceleration conclusions or model
adequacy. Seven fits have lower FreeRate likelihoods, by 0.1805–17.3648 log
units. They remain in the output and are listed in
`metadata/esmfold_freerate_optimization_flags.tsv`. Because the more flexible
model should be able to represent the Gamma categories, these are optimization
diagnostics requiring investigation before likelihood-based interpretation.
No flagged fit was silently replaced or removed. No model selection or
acceleration test is claimed. Full summaries, receipts and readback are tracked
under `metadata/esmfold_freerate_*` and `metadata/esmfold_rate_heterogeneity_*`.
