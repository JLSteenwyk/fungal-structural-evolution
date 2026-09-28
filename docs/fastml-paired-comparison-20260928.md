# Full paired FastML numerical comparison

Both numerical variants have completed all 156 coded-indel inputs, including
three empty inputs each. The full comparison verifies 9,492 source artifact
hashes and compares all 8,058,340 nonempty paired marginal entries at matching
node-descendant identities and character positions. Raw probabilities remain
unclipped, including the previously documented binary64 boundary excursions.

The cache-refresh variant changes fitted parameters in eight inputs, all in
OG0000294 (the exploratory Jaapia case). In total, 294 native probability
entries change; the maximum absolute change is 0.113505. No entry crosses a
presence probability of 0.5. This threshold comparison is descriptive and does
not imply stable uncertainty or a biologically validated insertion/deletion.

Two alignment inputs have worse independently replayed likelihoods after
cache refresh (about 0.186856 log-likelihood units); other affected fits improve.
Refreshing this cache does not ensure a better optimum. Preserve both variants
and require optimizer checks before selecting a result.

Using the original independent SciPy-rate replay, the largest absolute
reported/replayed likelihood discrepancy falls from 0.473948 to 0.003452 after
cache refresh. Seventy-four refreshed fits still exceed a descriptive absolute
difference of 1e-5. This threshold is not a fit-acceptance criterion.

A separate native-rate export links a small diagnostic program to the exact
frozen `libEvolTree.a` archive shared by both FastML variants. Across all 306
nonempty fits, its four-category gamma rates differ from the SciPy construction
by more than 1e-8, with a maximum category difference of 6.9184e-5. Mean rates
remain one within 2.23e-16. Native gamma discretization is therefore a candidate
explanation for some remaining replay discrepancies; rate comparison alone
does not prove that explanation.

The native-rate pruning replay has now completed all 306 fits. Its separate
audit verified 10,526 source, input, configuration and output files, the exact
fit set, exported rate identities and scalar arithmetic. For the 153 cache-
refreshed fits, the maximum likelihood discrepancy drops from 0.00345169 to
1.81917e-8, and the maximum marginal probability discrepancy drops from
6.80132e-5 to 1.78239e-9. Thus gamma discretization accounts for the previously
remaining discrepancies to this observed numerical precision.

The precision-only variant still has exactly three likelihood discrepancies
larger than the descriptive 1e-7 threshold: the previously identified OG0000294
envelope inputs, with errors of 0.469813–0.473949. This separates the stale-cache
problem from differences in rate calculation. Independent replay with native
rates reproduces the implemented model; it does not show that its gamma
approximation is mathematically exact or biologically preferable.

The two worse alignment optima remain worse under native-rate replay by
0.1868566884 log-likelihood units. Numerical reproduction therefore leaves an
optimization problem to resolve. No fit is automatically adopted from this audit.

The replay completed with one CPU core, about 106.5 CPU seconds and a 644 MiB
memory peak, without GPU or paid resources. A reproducible two-panel figure
shows likelihood and marginal-probability errors before and after matching
rate discretization. Its PDF and PNG are in
`results/figures/fastml-numerical-replay-20260928-v1`. Plot values below 1e-16
are displayed at that floor only; analytical values remain unmodified.

Reproducibility records:

- `scripts/compare_full_fastml_variants.py`
- `metadata/full_fastml_paired_comparison_20260928.json`
- `results/ancestral/full-fastml-paired-comparison-20260928-v1/paired_fit_comparison.tsv`
- `scripts/export_fastml_gamma_rates.cpp`
- `scripts/check_fastml_native_gamma_rates.py`
- `metadata/fastml_native_gamma_rate_checks_20260928.json`
- `metadata/fastml_native_rate_replay_plan_20260928.json`
- `metadata/fastml_native_rate_replay_launch_20260928.json`
- `scripts/audit_fastml_native_rate_replay.py`
- `metadata/fastml_native_rate_replay_audit_20260928.json`
- `scripts/plot_fastml_numerical_replay.py`
- `metadata/fastml_numerical_replay_figure_20260928.json`

Optimization, missing-data ascertainment and joint sequence/indel uncertainty
remain unresolved. None of these numerical comparisons qualifies an ancestral
structural ensemble or establishes a biological event.

## Effective optimizer limits and refinement preparation

A full audit of both native variants found that **all 306 nonempty fits** used
optimization level `low`, with effective limits of one outer cycle and one
model iteration. All 306 logs contain the model iteration-limit message.
Source inspection confirms that `updateOptimizationLevel(low)` overrides
both limits to one. Thus the optimization concern extends beyond the two
fits whose attained likelihood worsened. Successful numerical replay verifies
calculation at the returned parameters; it does not verify adequate fitting.

All 153 paired fitted trees are byte-identical across variants. A complete
refinement input set now preserves all 156 designs and supplies five starts
per input: both previous fitted parameter vectors and three fixed
shape/gain–loss-ratio settings. This yields 765 nonempty fits and 15 explicit
empty dispositions. Existing fitted trees are copied without rounding;
initial branch rescaling and empirical overwriting of the parameter starts
are disabled. Stationary mean-rate normalization is preserved.

Prepared controls use `mid` (which leaves explicit settings intact), 100
outer/model iterations and 1e-6 model/cycle tolerances. Effective native settings
must be checked at execution, alongside tree identity, starting values,
termination messages, candidate likelihoods and agreement across starts.
Two previous parameter starts can coincide; they are retained with provenance
and are not described as independent starts when identical.

The preparation is **not launched** at this checkpoint. Its resource plan
proposes eight CPU workers, 48 GiB RAM, no swap, 4 GiB per native process,
128 GiB output allowance and 12-hour per-fit timeouts. The broad 12–336 hour
processing allowance is uncalibrated because previous runs had only one model
iteration. No GPU or paid resources are planned.

Evidence:

- `scripts/audit_fastml_optimizer_settings.py`
- `metadata/fastml_optimizer_settings_audit_20260928.json`
- `scripts/prepare_fastml_optimizer_refinement.py`
- `metadata/fastml_optimizer_refinement_preparation_20260928.json`
- `results/ancestral/fastml-optimizer-refinement-inputs-20260928-v1/jobs.json`

The next step is a restartable full-grid runner with native effective-option
checks, independent numerical replay and explicit failed/unfinished outcomes.
No ancestral estimate is promoted by this preparation.
