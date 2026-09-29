# Complete unconstrained codon refits from profile-derived starts

The new run covers every one of the 1,632 local MG94 cases and all eight saved
profile solutions per case: 13,056 unconstrained optimizations. It follows the
complete profile audit, which found ten cases where constrained grid points
outperformed the single unconstrained reoptimization. Cases without an initial
flag receive the same full eight-start treatment.

Each start loads the original unconstrained model, then assigns all saved
branch and nuisance parameter values. The fixed target-branch constraint from
the profile export is not imported. Before optimization, a direct likelihood
calculation must reproduce the audited starting likelihood within 1e-6. The
run uses optimization precision 1e-7, with the same executable and fixed random
seed as the prior work. The final likelihood must not fall below its starting
value by more than 1e-5.

Every fitted model is exported with free parameter assignments and unchanged
non-parameter model content. A separate fresh HyPhy process reloads it and
checks its likelihood within 1e-6. Per-case receipts preserve all eight results,
starting and final likelihoods, parameter changes, target parameter, omega,
between-start deficits and raw artifact hashes. Original profile outputs remain
unchanged. A better result is a numerical improvement, not selection evidence.

The separately queued audit waits for the exact producer identity and successful
terminal state. It checks all source/seed/model identities, assignments and free
constraints; all saved parameter values against native logs; initial, optimized
and fresh-readback likelihoods; every aggregate table row; and all eight starts
per case. Its case table retains the spread between final starts and gains over
the prior unconstrained and best evaluated fits. Persistent differences remain
optimization concerns rather than being hidden by selecting a single result.

Producer: `scripts/run_local_mg94_multistarts.py`.
Auditor: `scripts/audit_local_mg94_multistarts.py`.
Plan: `metadata/local_mg94_unconstrained_multistart_plan_20260927.json`.
Launches: `metadata/local_mg94_unconstrained_multistart_launch_20260927.json` and
`metadata/local_mg94_unconstrained_multistart_audit_launch_20260927.json`.
Outputs: `results/cds/local-mg94-unconstrained-multistarts-20260927-v1` and
`results/cds/local-mg94-unconstrained-multistart-audit-20260927-v1`.

The producer has four one-CPU workers, an 8 GiB memory limit, no swap, and a
20 GiB output allowance. The planning range is 1–48 hours, based on the previous
10 CPU-hour profile batch with additional uncertainty for tighter unconstrained
optimization. The audit uses one CPU, 8 GiB memory and 0.1–4 hours after fitting.
These are planning allowances, not measured completion forecasts. Startup is
confirmed by advancing completed-case checkpoints and four live HyPhy workers.

The full results and audit are pending. Multiple starts do not establish a
global optimum, synonymous-distance confidence intervals, absence of saturation
or eligibility for selection inference. Biological interpretation and selection
tests remain separate stages. GPU prediction remains paused.


## Partial completed-case disagreement review

The frozen September 27 snapshot retains all 1,632 cases: 777 complete and
855 pending. All eight starts for each completed case were retained, their
saved artifact hashes checked, and starting/optimized/fresh-readback score
arithmetic verified. This does not replace the queued full model-constraint,
parameter and source-identity audit.

Among completed cases, between-start log-likelihood spread exceeds 0.00001
in 90 cases, 0.01 in 75, 1 in 22, and 10 in 14. The largest spread is
162.658494 for `Malassezia__5000823at2759__code1`: seven starts agree within
0.00001 of the best likelihood (omega approximately 0.079975), while the
factor-10 start remains much poorer (omega approximately 0.004166).
The next largest spreads are 120.768979 for marker 4813267at2759 and
77.145974 for marker 4987397at2759, also in Malassezia/code1 and with seven
near-best starts. These examples show optimizer sensitivity rather than
biological evidence for selection. They do not establish global optima.

No completed case has a near-best-start omega range greater than 0.1 under
the stated 0.00001 likelihood tolerance. This coarse descriptive threshold
is not a confidence interval or proof of identifiability. Partial completion
can be biased toward faster cases, so these counts are not full-grid rates.
Summary counts and maxima were independently recomputed from the frozen case
table. Evidence: `metadata/local_mg94_multistart_partial_review_20260927.json`;
full case dispositions remain outside Git in
`results/cds/local-mg94-multistart-partial-review-20260927-v1`.

## Full parameter-range summary queued September28

After the complete1,632-case/13,056-fit audit reaches terminal success,
`summarize_mg94_multistart_parameters.py` will retain the minimum and maximum
of every saved free parameter across all eight starts and across starts
within1e-5 log-likelihood of the best observed start. This uses the audit's
existing near-best tolerance; it does not choose a new biological eligibility
threshold. Branch and nuisance parameters remain explicit alongside omega.
Near-best parameter spread is descriptive, not a confidence interval or proof
of identifiability, saturation, a global optimum or selection.

Every case/fit checksum and source-audit binding is checked. Counts of near-best
starts and best likelihoods must agree with the complete audit. All emitted
range arithmetic and enclosing bounds are read back. The reproducible fixture
`check_mg94_parameter_ranges.py` retains parameter variation at near-equal
likelihood, preserves all-start extrema, excludes poorer fits only from the
near-best range, handles one best start and rejects nonfinite values.

Plan `metadata/mg94_parameter_range_plan_20260928.json` and launch record
`metadata/mg94_parameter_range_launch_20260928.json` bind the inputs and live
service. Resources are one CPU,8GiB RAM, no swap,1GiB output allowance and
0.1–4 active hours excluding wait. No new optimization or GPU prediction is
launched. Output will be
`results/cds/local-mg94-multistart-parameter-ranges-20260928-v1`.
The source audit and this summary remain pending at launch.

## Full run and audit completed September 28

All 1,632 cases and 13,056 unconstrained fits finished, followed by successful
full artifact/numerical audit. The audit checked 65,280 artifact hashes and
224,752 parameter values; the maximum fresh saved-likelihood discrepancy was
2.365e-11. Completion bindings and all three successful terminal states are in
`metadata/local_mg94_multistarts_completed_20260928.json`.

Between-start likelihood spread exceeds 1e-5 in 204 cases. The final long-running
Malassezia case (`Malassezia__649304at2759__code1`) has a spread of 0.033685,
three starts within 1e-5 of its best observed log likelihood, and an improvement
of 0.015709 over the prior unconstrained result. This supports retaining
optimization uncertainty; it does not establish a global optimum.

The full parameter-range stage also finished, retaining 28,094 parameter rows
across all cases, with ranges over all starts and over starts within 1e-5 of the
best observed likelihood. Source-file bindings and serialized range arithmetic
were checked by production. An independent reconstruction of those ranges
remains outstanding. These numerical ranges are not confidence intervals or
selection evidence, and do not resolve saturation or parameter identifiability.
