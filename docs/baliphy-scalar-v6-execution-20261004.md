# Complete corrected scalar sampler execution

The full corrected startup grid has closed successfully: 1,620 native roles,
405 four-chain groups, 135 effective inputs and 324 original configuration
aliases. Every startup passed. Independent serialized readback and both
original process journals are included in the 20,397-binding completion archive.
The [startup closure](../metadata/baliphy_scalar_v6_preflight_completed_20261004_v1.json)
does not establish MCMC convergence or ancestral uncertainty.

The corrected native short sampler and concurrent read-only resource observer
are now running. The [execution plan](../metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json)
retains every role, all three original priors and the existing 20-iteration
computational horizon. Its copied jobs are byte-identical to the complete
software-qualified future grid. No historical numeric values are adopted.
The [preparation receipt](../metadata/baliphy_scalar_v6_execution_preparation_20261004_v1.json)
and [original terminal proof](../metadata/baliphy_scalar_v6_execution_preparation_transport_20261004_v1.json)
bind the actual preparation. Original wait 54212 exited zero; the exact wrapper
identity, whole initial/terminal journal payloads and manager start/end match.

Six original controllers handle sampler production, independent sampler replay,
sampler closure, concurrent telemetry, independent telemetry replay and telemetry
closure. The [launch inventory](../metadata/baliphy_scalar_v6_sampler_execution_launches_20261004_v1.json)
records them. Admission waited for the original corrected startup closure and
all four earlier historical prerequisite closures. The observer waits for
startup only and watches the sampler while it runs; it does not depend on
sampler completion. Its own actual receipt and closure belong after native
termination. There are no automatic retries or replacements of failed jobs.

The sampler has 16 workers and 16 CPUs, 200 GiB cgroup memory, no swap and one
BLAS thread per native process. Declared address-space leases total at most
192 GiB, leaving 8 GiB cgroup headroom. There are 1,512 roles with 12 GiB leases
and 108 with 48 GiB leases. The reservation budget can reduce simultaneous
workers for the larger roles. Reader/closer scopes have two CPUs and 32 GiB;
the observer has one CPU and 2 GiB, and its reader/closer have two CPUs and 4 GiB.
Available RAM and disk were freshly checked at launch. The resource plan's
inherited historical RAM availability is not a current machine reading.

Output planning allows 128 GiB with a 256 GiB free-disk admission guard;
per-file/native address-space/time limits remain unchanged. The sum of role
wall timeouts is a bound on declared worker budgets, not a finish ETA.
Actual telemetry records original PID/create/command identities, native limits,
process memory readings and cgroup counters. Missed, fast, transitioning and
failed attempts remain explicit. Cgroup memory includes descendants/cache and
is not a native RSS peak or guaranteed exact per-role maximum.

At 08:55 UTC, the [original runtime checkpoint](../metadata/baliphy_scalar_v6_sampler_execution_checkpoint_20261004_0855.json)
records 16 completed dispositions: 15 finite integrity checks and one explicit
special-value review, with 20 original attempt identities. Neither full
sampler nor observer has a producer receipt, readback or final closure yet.
Early counts are unclosed progress, not completed full-grid analysis.

The reviewed role is OG0000230, broad prior, chain 3. At iteration 1 its native
scalar JSON explicitly tags positive infinity for `ASRV.Gamma:alpha`, and the
paired TSV prints `Infinity`; other iterations include very large finite alpha.
The role's exit is zero, but it is not numerically admitted. Original logs and
all source inputs remain unchanged, with no exported ancestral arrays. The
special-value validator separates this observed state from malformed JSON and
from the installed formatter's exponent bug. Its proposal/model cause is not
established; preserve the complete grid before proposing model changes.
The [independent completed-role reconstruction](../metadata/baliphy_scalar_v6_special_value_role_evidence_20261004_v2.json)
matches the entire saved disposition and checks all 21 native scalar rows,
retaining the exact JSON and TSV special-value tokens. A new V1 report reader
initially used a normalized section name as a serialized MCON key and stopped
with `KeyError`; its source and failure receipt remain. The separate V2 report
reader uses the actual `parameters//` field. This reporting failure did not
affect the production validator, original job or native outputs.

Reproduction uses the preparer and launcher:

```text
scripts/prepare_baliphy_scalar_v6_sampler_execution_v1.py
scripts/launch_baliphy_scalar_v6_sampler_execution_v1.py
scripts/record_scalar_v6_sampler_execution_checkpoint_v1.py
scripts/record_scalar_v6_special_value_role_evidence_v2.py
```

Launched plans, source bytes and output roots are immutable. Do not rerun the
preparer/launcher into these namespaces. A new execution requires a new version,
fresh admission evidence and explicit preservation of all earlier outcomes.

Full reader/source/artifact/original-journal closure remains mandatory. Twenty
iterations are insufficient for accepted posterior uncertainty. Long chains,
mixing/ESS/convergence, prior/model adequacy, predictor/root controls and the
eight evolutionary aims remain open. GPU prediction stays paused; no paid
resources have been provisioned.
