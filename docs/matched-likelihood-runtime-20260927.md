# Likelihood runtime measurements before the full fitting grid

The current likelihood evaluator has been timed on synthetic numerical inputs
whose record/background/family counts come from the complete verified record
summaries. It uses the full 242-column species factor and five fixed-effect
columns, with all numerical-library thread pools restricted to one. This is a
computational benchmark, not a biological pilot or a fit to selected project
results.

The record-count median across settings is 1,650, but there is no setting of
that size: its nearest bracketing observations are 980 and 2,320. Both are
measured, along with the minimum and maximum. Each size receives five timed
repeats under each of three variance-component settings, for 60 measurements.
Each timing includes fresh covariance construction and a profiled restricted
likelihood calculation; data generation and warm-up are outside the timings.

| Records | Full phylogenetic evaluation, median wall time |
|---|---|
| 148 | 2.9 ms |
| 980 | 14 ms |
| 2,320 | 38 ms |
| 10,963 | 160 ms |

These measurements reflect the current host while other authorized CPU jobs
are active. Raw wall/CPU times, ranges, library/thread details and source hashes
are in `metadata/matched_likelihood_runtime_20260927.json`. The script is
`scripts/benchmark_matched_likelihood.py`.

The three previously checked synthetic optimization fixtures required 367,
768 and 790 objective calls across their 22 boundary/start candidates, plus
seven final readback/gradient calls each. Those counts are not a bound on
convergence costs for real data. Some evaluations omit the species term and
are much faster. Full-grid elapsed time also depends on exact-input reuse,
parallelism, per-setting rank, conditioning, I/O and optimizer review cases.
No whole-project or full-grid completion ETA is inferred from these timings.

At 414,720 fits before exact reuse, repeatedly recomputing invariant group
summaries and matrix products can consume substantial CPU time. The next
implementation step is a cached evaluator that precomputes unchanged sufficient
statistics, while retaining the checked evaluator as its reference. Numerical
equivalence across boundaries and cancellation-sensitive cases, followed by
repeat timing, is required before adopting it for production fitting.
