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
