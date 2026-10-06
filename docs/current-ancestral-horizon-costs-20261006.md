# Current full ancestral run costs and longer-horizon requirements

The current resource census and its independent full readback have completed
across **1,620 original roles, 135 effective inputs and 405 prior quartets**.
All **24 original native failures** and their **24 separate V10 diagnostic
comparisons** remain explicit: 1,644 measured attempts in total. Every one of
**16,440 native files / 2,784,481,010 bytes** is freshly hashed, sized and checked
against its closed original source. The independent reader reconstructs all
file measurements, attempt dispositions and cost arithmetic.

This updates the older 1,000-iteration cost reference to the current logger
and native outcomes. It does not run a pilot or reduce the 501-fungal-entry
plus 25-outgroup study. These focal ancestral configurations support one part
of the full project; their protein tips, aliases and chains are not independent
biological samples. No new native sampler or GPU inference runs here.

## Observations and the computational cost baseline

All 1,620 original V6 attempts and all 24 V10 comparisons are measured.
All measured attempts consumed **37.28 aggregate worker-wall hours**, including
the separate comparisons. For a prospective cost scenario, the 1,596 original
native-zero outcomes and the 24 native-zero V10 comparisons supply exactly
one computational baseline per original role. The original failed attempts
remain in the attempt and file ledgers; they are not silently replaced in
the scientific dataset. These cost baselines total **36.49 worker-wall hours**.
Worker wall time is not a direct CPU-time measurement.

Fourteen baselines retain special-value review dispositions. They are included
for resource accounting, with scientific and posterior eligibility explicitly
false. Memory observations remain tied to their original V6 processes and
sampling gaps. Per-native live V10 memory was not observed, and short-run memory
does not establish a safe long-chain peak.

The [complete public 1,644-attempt table](tables/current_ancestral_horizon_attempts_20261006_v1.tsv)
is a byte-identical copy of the verified source table. The detailed 16,440-file
ledger remains outside Git at
`results/ancestral/current-full-horizon-costs-20261006-v1/native_files.tsv`;
its hash and all input locations are recorded in the
[producer receipt](../metadata/current_ancestral_horizon_cost_20261006_v1.json).

## Explicit 10,000-iteration scenarios

The existing proposed first long horizon is 10,000 iterations. Native scalar
and latent traces would contain 10,001 saved states; with the unchanged
ten-iteration frame interval, alignments and joint site-property outputs would
contain 1,001 saved frames. Their multipliers differ: **10,001/21 scalar states**
and **1,001/3 frames**. Fixed artifacts keep their measured byte counts in the
scenario. Every byte projection uses exact integer ceiling arithmetic.

The resulting existing-native-file projection is **935,896,642,895 bytes**.
An additional **9,577,915,200-byte** latent-output allowance uses the largest
observed 600-byte physical diagnostic line for the 1,596 baselines without that
new file. Together these give **945,474,558,095 bytes (0.945 TB)** of native
outputs under the stated assumptions. This omits derived projection arrays,
later summaries, temporary/readback products and potential alignment growth;
it is not the complete disk requirement or a proven output bound.

Scaling the current worker-wall measurements by 10,000/20 gives approximately
**18,246 aggregate worker-wall hours**. Dividing that total by workers gives:

| Workers | Ideal allocation of the linear scenario |
| ---: | ---: |
| 4 | 190.1 days |
| 8 | 95.0 days |
| 16 | 47.5 days |

These are conditional planning calculations, **not finish or convergence
estimates**. Initialization overhead, changing alignments, effective worker
availability, memory leases, logging and later validation can change runtime
in either direction. The older successful-only 1,000-iteration reference scales
to about 16,121 worker hours, with two failed roles unestimated; that older
implementation does not independently calibrate the new horizon.

## File-cap implication and remaining launch requirements

**96 projected files exceed the existing 2 GiB per-file cap**. The largest
linear file projection is 12,955,334,726 bytes (12.07 GiB), so a prospective
16 GiB cap covers this scenario while retaining some margin. It does not bound
longer-chain growth. The
[complete file-cap assessment](../metadata/current_ancestral_horizon_output_budget_assessment_20261006_v1.json)
records that calculation; it changes no live configuration or limit.

Before a long native horizon is admitted, the workflow needs full-horizon
construction/readers, a census of derived-output costs, explicit memory/CPU/
file/disk bounds and aggregate disk monitoring. Original failures and numeric
reviews must remain explicit. Fresh chains need their own seeds and output
namespace; this census does not resume, concatenate or accept an earlier chain.
No long-horizon sampler is launched or queued by this stage. Convergence,
model adequacy, accepted genealogies, structural prediction uncertainty and
all eight evolutionary aims remain unfinished.

## Execution and reproduction

Original producer session **87251**, invocation
`d4d565feceb543f3956a8dece6ed77a2`, exited zero in 27 seconds. Original independent
reader session **13397**, invocation `74b173b2f191488c911d7403e52df52d`, exited
zero in eight seconds. Complete original wrapper messages, invocation start/end,
used-source and output bindings close separately for both:

- [Producer execution closure](../metadata/current_ancestral_horizon_cost_transport_20261006_v1.json).
- [Independent full result](../metadata/current_ancestral_horizon_cost_readback_20261006_v1.json).
- [Reader execution closure](../metadata/current_ancestral_horizon_cost_readback_transport_20261006_v1.json).

Both used two CPU equivalents, a 32 GiB/no-swap cgroup, 24 GiB address-space
limit and one BLAS thread. Two-hour CPU and four-hour wall bounds are safety
caps. The eight-GiB output allowance is a planning reserve; a 256 MiB per-file
limit is enforced. No paid resources are used. Original native inputs, existing
closure records and SHA helpers are shared dependencies; the reader uses a
separate chunked counter, divmod byte arithmetic and Decimal cost reconstruction.

The producer command was:

```bash
python scripts/inventory_current_ancestral_horizon_costs_v1.py \
  --output results/ancestral/current-full-horizon-costs-20261006-v1 \
  --receipt metadata/current_ancestral_horizon_cost_20261006_v1.json \
  --iterations 10000
```

The independent command was:

```bash
python scripts/readback_current_ancestral_horizon_costs_v1.py \
  --producer metadata/current_ancestral_horizon_cost_20261006_v1.json \
  --producer-transport metadata/current_ancestral_horizon_cost_transport_20261006_v1.json \
  --receipt metadata/current_ancestral_horizon_cost_readback_20261006_v1.json
```

Reproduction needs the original raw native files and source closures, the
recorded bounded wrapper/cgroup and fresh receipt/output namespaces. Scripts
and resource plans are versioned; the large native corpus and file ledger
remain outside Git. Raw-data public release remains outstanding.
