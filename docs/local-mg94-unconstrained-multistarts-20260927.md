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
