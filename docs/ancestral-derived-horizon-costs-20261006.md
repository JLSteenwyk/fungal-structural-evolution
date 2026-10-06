# Full ancestral derived storage census

The full current derived-output census and independent readback are complete
across **1,644 attempts**: 1,620 original V6 roles and 24 separate V10 diagnostic
comparisons. All 24 original native failures and 14 numerically reviewed cost
baselines remain explicit. This stage creates no new sampled or predicted
structures and does not establish an adequate ancestral posterior.

## Measured current outputs

All **6,477 non-native files / 1,813,085,268 bytes** in the two completed run
trees were freshly hashed and sized. Native attempt directories are excluded
because they have their [separate complete census](current-ancestral-horizon-costs-20261006.md).
All 4,818 exported projection files contain exactly ten expected NumPy members.
The producer streams NPY member payloads; the independent reader loads NumPy
values and checks every shape, dtype and raw-data digest against the original
closed frame records. All **48,180 members** pass. The reader independently
measures NPY headers from their version/length prefix and reconstructs each
file's serialization overhead. Both paths use NumPy serialization as a shared
dependency; this does not repeat the earlier native joint-frame/likelihood checks.

| Current derived product | Files/bytes |
| --- | ---: |
| Exported projection files | 4,818 / 1,432,750,716 bytes |
| Raw array payloads | 1,419,886,656 bytes |
| NPY/ZIP serialization overhead | 12,864,060 bytes |
| Per-role disposition files | 100,159,817 bytes |
| Stage metadata, including proof archives | 280,174,735 bytes |
| Largest individual projection file | 5,702,962 bytes |

The public [complete 1,644-attempt table](tables/ancestral_derived_horizon_attempts_20261006_v1.tsv)
is a byte-identical copy of the verified source. Detailed projection-file,
array-member and complete non-native-file ledgers remain outside Git under
`results/ancestral/full-derived-horizon-costs-20261006-v1/`. Their paths and
hashes are in the [producer receipt](../metadata/ancestral_derived_horizon_cost_20261006_v1.json).
Zero-byte lock files are measured, but are not evidence that a process or
analysis completed. The V10 supplementary root receipt is freshly measured;
scientific custody comes from the original closed disposition/frame sources.

## Explicit 10,000-iteration disk scenario

The proposed frame interval remains ten iterations: 1,001 saved frames per
role. For all 1,606 cost baselines with integrity-checked exports, each role's
three observed file sizes is scaled by 1,001/3 using exact integer ceiling
arithmetic. That gives **478,061,155,643 bytes** of array files.

The remaining 14 numeric-review baselines have no admitted exports. For resource
planning only, each uses the largest measured single frame from an
integrity-checked peer with the **same effective input and prior**, multiplied
by 1,001. Donor identity and size are explicit in the full table. These labelled
proxies add **16,764,325,578 bytes**. They are not reconstructed draws,
biological imputations or evidence that the original reviewed arrays are valid.
All original failed attempts stay in the ledger, with zero scenario contribution;
their separate V10 comparison contributes the single computational baseline.

| Planning component | Projected bytes |
| --- | ---: |
| Observed-role arrays | 478,061,155,643 |
| Labelled resource proxies for 14 missing-export roles | 16,764,325,578 |
| All arrays under these assumptions | 494,825,481,221 |
| Previously verified native outputs plus latent allowance | 945,474,558,095 |
| Native plus arrays | **1,440,300,039,316** |

This is approximately **1.44 TB in decimal units**, before long-chain metadata
growth, temporary/readback products, posterior summaries and filesystem
allocation overhead. Current metadata bytes are measured but are not projected
into longer runs. Alignment changes and extreme native states can change future
file sizes. Neither the observed-role scenario nor the peer maximum is a proven
upper bound. No live file cap or process configuration is changed here.

Full-horizon construction and streaming readers, long-chain memory/CPU/file
limits and aggregate disk monitoring remain required before a long sampler can
launch. Mixing, convergence, model adequacy and accepted biological frameworks
remain separate requirements. All eight evolutionary aims are unfinished.

## Original execution, failure retention and reproduction

The bounded producer's original API session **61167** and successful V2 reader's
original session **27623** both exit zero. Their exact original wrapper messages,
invocation start/end and all used-source/output bindings close independently:

- [Producer original execution closure](../metadata/ancestral_derived_horizon_cost_transport_20261006_v1.json).
- [Full independent result](../metadata/ancestral_derived_horizon_cost_readback_20261006_v2.json).
- [Reader original execution closure](../metadata/ancestral_derived_horizon_cost_readback_transport_20261006_v2.json).
- [Completed handoff](../metadata/ancestral_derived_horizon_cost_completed_20261006_v1.json).

The initial reader session **17168** exited one before array checking because
it expected a direct completion-archive label. The producer had correctly found
that archive already closed in a preceding source transport. The exact failed
source, payload, stderr and invocation journal are
[preserved](../metadata/ancestral_derived_horizon_cost_readback_failure_20261006_v1.json).
A fresh V2 reader recognizes either closed lineage while still requiring the
exact completion-archive hash. That is the only source change; numeric,
array-value, arithmetic and eligibility checks are unchanged. The producer and
native samplers were not restarted.

Both successful stages used two CPU equivalents, a 16 GiB/no-swap cgroup,
12 GiB address-space ceiling, one BLAS thread, two-hour CPU/four-hour wall
safety limits and a 256 MiB per-file limit. The eight-GiB output allowance is
a planning reserve. Actual elapsed times were about 31 and 14 seconds;
they measure census work, not the duration of future MCMC runs. GPU prediction
remains paused and no paid infrastructure is used.

```bash
python scripts/inventory_ancestral_derived_horizon_costs_v1.py \
  --output results/ancestral/full-derived-horizon-costs-20261006-v1 \
  --receipt metadata/ancestral_derived_horizon_cost_20261006_v1.json

python scripts/readback_ancestral_derived_horizon_costs_v2.py \
  --producer metadata/ancestral_derived_horizon_cost_20261006_v1.json \
  --producer-transport metadata/ancestral_derived_horizon_cost_transport_20261006_v1.json \
  --receipt metadata/ancestral_derived_horizon_cost_readback_20261006_v2.json
```

Existing result namespaces are immutable; reproduction requires fresh output,
receipt, execution and unit names. Full dataset dimensions, original failures,
source hashes and the stated scenario must remain unchanged.
