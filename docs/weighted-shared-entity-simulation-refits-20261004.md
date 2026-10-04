# Dependent Gaussian responses with actual frozen refits

The response generator now connects to the actual working-model optimizer and
independent numerical reader. This addresses a missing step between dependent
simulation and future uncertainty calibration. It does not replace full source
admission, complete numerical/timing/fitting closures, biological model checks,
or the eight evolutionary aims.

## Source and model contracts

`scripts/refit_weighted_shared_entity_simulation.py` requires the original fit
draft SHA, original upstream candidate identity SHA, original response SHA and
the upstream Gaussian model contract. The latter binds the ordered selected
labels, canonical signed sparse incidence operators, species factor F,
positive residual diagonal D, design X and exact component names. This prevents
a change in coefficient interpretation from slipping through a covariance
Gram check that only depends on the column space of X. Original coefficient
mapping lengths are checked when supplied by the source. The production
whole-grid source reader must still reconstruct the actual coefficient and
row mapping; this per-model bridge cannot confer source closure.

The original identity, covariance audit, route, response and upstream contracts
stay distinct from the new simulated response identity. Only the selected
cohort's arrays are copied; rebasing its row indices avoids copying the global
factor bank for each candidate. Ordered original row indices remain hashed.
The immutable original fit draft retains all four reuse policies, signed and
unsigned modes, five trees, two methods and both outcomes. Production settings
are unchanged: three original optimizer starts, 500 evaluations, 200 iterations,
gradient tolerance 1e-6 and maximum scaled variance 1e6; independent replay
objective tolerance is 1e-7 and gradient tolerance is 1e-6, with column batch 32.
Generating variance ratios must lie inside that same scaled box. Generating
truth never becomes an optimizer starting point.

The calibration contract includes the exact source/model identity, immutable
fit draft SHA, optimizer and independent-reader settings, interval definition,
code hashes and interpreter/package/platform versions. Existing qualified
parents and math scripts are not changed. Exact response replay is required
in the recorded environment; portability to a different build is not asserted.

## Nominal interval target and unresolved outcomes

The prespecified target is the normal Wald interval with nominal probability
0.95, fitted conditional coefficient covariance and the standard-normal
0.975 quantile. The same rule applies to ML and REML. SciPy's
[normal distribution and inverse CDF documentation](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.norm.html)
defines the quantile calculation. Choosing this conditional working interval
is a project calibration target, not a claim that estimated-variance uncertainty
is covered or that its biological coverage equals 95%. The frozen methods
retain their own ML/REML scale calculations; no post-fit reference change is
used to improve apparent coverage.

Only `independently_audited_working_candidate_pending_inferential_calibration`
from the frozen reader can yield
`refit_independently_checked_pending_calibration`. Other numerical outcomes
remain `refit_requires_review` or `refit_error_requires_review`. Reviewed or
failed attempts stay in every prescribed replicate denominator. The frozen
reader reproduces primary numerical failures rather than accepting invented
failure records. Invariant assertions remain fatal and preserve the original
attempt inputs and primary payload; they are not converted to soft successes.

Saved refit replay regenerates the complete Gaussian stream/response record,
reconstructs the simulated identity, runs the frozen independent reader on the
serialized primary result and recalculates the interval/status fields. Coverage
accounting first replays all records and refuses foreign source, method,
settings or interval contracts. Every prescribed replicate is required. Its
exact fixed-sample Monte Carlo envelope retains unresolved outcomes; this is
not a coefficient interval, simultaneous coverage statement or acceptance of
the biological model.

## Verification design and original attempt

`scripts/check_weighted_shared_entity_simulation_refits.py` declares 128 actual
unmocked refits: the 64 original synthetic model axes crossed with zero and
positive generating shared-effect variance. The axes cover both loading modes,
all four reuse policies, both outcomes, ML/REML and pair-exception cases; q4,
q5 and q6 covariance bases are included. Each primary result is saved before
readback. Generated responses, source identities, primary payloads, independent
checks and accounting records are stored outside Git in
`data/software_audits/weighted-shared-entity-refits-20261004-v1/`.

Every serialized record is independently replayed, and selected fitted points
are compared with direct dense Gaussian likelihood/coefficient/covariance/
scale calculations. Private-copy source/model, seed/truth, identity, optimizer
history, readback, interval and accounting alterations must be rejected. A
separate explicit injected-reader-error fixture checks the unresolved branch;
it is never included in the count of actual numerical refits. Published parent
files are never temporarily changed for corruption tests.

Resources are two CPUs, 16 GiB, no swap, one BLAS thread, 12 GiB address space,
3,600 CPU seconds, 7,200 wall seconds and 256 MiB per file. Minimum free RAM/disk
are 16/128 GiB; output allowance is 2 GiB, with no GPU or new charge. The original
wrapper is managed by `systemd-run --user --wait --collect`; actual tool-wait
exit, wrapper PID/create/command, invocation-specific manager start/end and
the entire original terminal payload must match before qualification.

Verification receipts use the prefix
`metadata/weighted_shared_entity_refit_software_` and suffix
`20261004_v1.json`: resources, execution, validation and transport. The
transport reader verifies all fixture/code/parent hashes and the original
terminal receipt SHA. Large response and numerical payloads remain outside Git;
the repository stores their hashes and exact paths.

The original software wait 3751 exited zero. All 128 primary fits had fitted
variance vectors and passed serialized numerical replay and dense selected-point
comparisons. Primary outcomes were 78 optimized candidates and 50 optimizer
reviews. Independent results retained 34 curvature/search reviews in addition
to the 50 producer reviews; 44 passed the independent working checks. All 84
reviews remain unresolved, and every record retains scientific eligibility
false. Maximum dense objective/coefficient absolute differences were
4.70e-10/3.77e-11. These numerical differences are not uncertainty estimates.

Verification rejected 15 copied-source/model alterations, 15 refit replay
alterations, four interval changes, two generating-ratio errors and three
accounting changes. The injected reader exception remained unresolved and
could not replay as an actual numerical result. One replicate per scenario
does not support any coverage estimate.

Original wrapper PID 123957, invocation
`416ce3bd97834bcaad6814163f95a775`, actual wait, exact initial/terminal messages
and manager start/completion were verified. Software CPU time was 124.34 seconds,
wall time 128.52 seconds and reported native child peak RSS 155,181,056 bytes.
Native RSS is distinct from cgroup memory accounting. The shared wrapper's
historical generic scope text is inherited; this validation and transport
describe the actual refit checker. No real production source was admitted.

## Remaining production work

One replicate per synthetic model/scenario verifies software integration. It
does not estimate coverage, implement the complete fungal calibration design,
or constitute a biological pilot. Complete the original full weighted
numerical and timing closures, reconstruct full source admission and run the
unchanged production fits with full independent output closure. Then specify
generating scenarios, fixed replicate counts, precision and compute/storage
estimates before biological simulation. Calibration must address model adequacy,
cross-model/global dependence, tree and ancestral uncertainty, selection and
multiple testing. The full structural atlas, accepted species/gene/reconciliation/
dating framework, all eight aims and experimental hypotheses remain required.
