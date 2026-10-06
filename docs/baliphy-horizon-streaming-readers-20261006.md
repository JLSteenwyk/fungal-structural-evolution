# Horizon-aware streaming ancestral readers

The new V2 reader separates a declared sampling horizon from the old
20-iteration qualification constants. Scalar rows and joint frames are read
incrementally, and export-reference ledgers can be iterators. Each frame is
decoded, checked and released before the next frame is processed. Memory still
depends on the largest raw frame, its decoded arrays and transient parser
allocations; this is not a qualified long-chain peak-memory bound.

The full software qualification replays **all 1,644 current native attempts**:
1,620 original V6 roles and 24 separate V10 comparisons. Original native
configuration/process/artifact custody and inherited outcomes are retained.
The 24 original native failures stay explicit. All 14 numeric reviews keep
their excluded arrays; they do not acquire new exports. The new code does
not change the qualified short-run modules, native models, seeds or priors.

## Trace and export contracts

- Scalar JSON headers, duplicate keys, explicit numeric-quality tags and
  iteration order are checked for every state from zero through the declared
  horizon. The primary and independent qualified scalar validators agree on
  each row's numeric census and review state. Reviews can be emitted one at
  a time through a sink rather than collected into a long in-memory list.
- Native TSV has the same complete horizon and unique column headers. Native
  prior, likelihood and posterior remain finite and satisfy the existing
  absolute `1e-7` sum check. Fully finite roles receive the original mapped
  field check at relative `2e-13` and zero absolute tolerance. Tagged
  nonfinite/literal-null roles remain reviews and receive no admitted arrays.
- Optional V10 latent traces are paired row-by-row with their original scalar
  records. The existing strict Decimal parser and 90-digit exponentiation,
  Laplace-density, category-rate and explicit overflow checks remain unchanged.
  Finite latent values explaining overflow do not remove scalar exclusions.
- Source and runtime trees retain their exact rooted clade correspondence and
  the unchanged strict `1e-10` branch-length check. Original source-node
  mappings, four candidate nodes, extant observations and negative-edge review
  counts are preserved.
- FASTA and joint traces require exactly the declared saved-frame sequence.
  The default interval is ten; a horizon that is not divisible by ten does
  not invent a final frame. Every joint frame passes the qualified all-node,
  sequence/category/rate/anchor decoder. Extant constraints and separate FASTA
  frame integrity are checked without asserting identity between separate
  logger draws.
- A readback checks every original exported array value, shape, dtype and
  source hash. Reference ledgers can be streamed, so longer runs do not require
  retaining every frame's reference metadata or decoded arrays together.

The full current corpus contains **34,020 scalar rows, 504 latent diagnostic
rows, eight explicitly retained overflow rows, 4,860 native FASTA frames and
4,818 previously exported joint frames**. Output integrity remains distinct
from posterior adequacy and biological framework acceptance.

## Transactional exports and software controls

`stage_exports` first checks scalar/latent traces and FASTA frames. A numeric
review creates no array namespace. Otherwise each frame is serialized and
verified in a fresh `.pending` directory with an incremental JSONL ledger.
Only after every frame passes and the raw input hashes remain unchanged does
the helper write its completion record and rename the pending directory.
An error retains the pending products and failure record, with no completed
export directory. Existing namespaces cannot be reused or overwritten.

`verify_staged_exports` independently reads the raw traces again, verifies
the sealed ledger/source hashes, checks exact file membership and compares
every serialized array against newly decoded raw frames. It streams the
reference ledger. Neither helper verifies native execution custody itself;
that remains a mandatory outer-controller requirement. Their scientific and
posterior eligibility flags remain false.

Three positive software controls cover a 17-iteration horizon, a
10,000-iteration horizon and an explicitly tagged numeric review. The long
control contains exactly **10,001 scalar states and 1,001 joint frames**, with
transactional export and full streamed readback. These are synthetic file
controls, not sampled ancestors, MCMC runs or a biological pilot.

Twenty-nine rejection controls include missing/extra/reordered scalar and
joint rows, duplicate JSON/TSV keys, unencoded nonfinite states, changed mapped
scalars, altered FASTA observations, latent density/prior errors, an invalid
late category, missing/foreign arrays and a self-consistently rehashed extra
ledger row. A late-frame failure retains one pending frame but no final
directory. The tagged review creates zero arrays.

## Execution and preserved failure

The initial original qualification wait **92824** exited one at source-binding
aggregation: it had not imported the existing original V10 readback transport
that closes the V10 attempt configuration. Its exact source, original payload,
stderr, wrapper identity and invocation journal are
[retained](../metadata/baliphy_horizon_streaming_reader_failure_20261006_v1.json).
The fresh V2 checker imports that existing transport. The V2 module additionally
streams export-reference ledgers, releases consumed frame objects and checks
sealed long-control arrays through their streamed ledger. Original numerical
tolerances and all current native data remain unchanged; no sampler is restarted.

Full qualification and original execution closure are recorded in the
[complete result](../metadata/baliphy_horizon_streaming_reader_qualification_20261006_v2.json),
[original transport](../metadata/baliphy_horizon_streaming_reader_transport_20261006_v2.json)
and [completed handoff](../metadata/baliphy_horizon_streaming_reader_completed_20261006_v2.json).
The public [full attempt readback table](tables/baliphy_horizon_streaming_reader_roles_20261006_v2.tsv)
preserves every original failure, review and separate comparison.

The qualification uses four CPU workers, a 32 GiB/no-swap cgroup, 24 GiB
address-space ceiling, one BLAS thread, a three-hour per-process CPU/six-hour
wall safety limit and a 256 MiB per-file limit. The eight-GiB output allowance
is a planning reserve. This reader allocation does not qualify future MCMC
memory, CPU, disk or file bounds.

```bash
python scripts/check_baliphy_horizon_streaming_readers_v2.py \
  --output data/software_audits/baliphy-horizon-streaming-readers-20261006-v2 \
  --receipt metadata/baliphy_horizon_streaming_reader_qualification_20261006_v2.json \
  --workers 4
```

Reproduction requires the pinned full source corpus and fresh result, receipt,
execution and unit names. Existing namespaces are immutable. These software
and full short-corpus checks do not validate a native 10,000-iteration horizon,
choose burn-in, prove convergence or admit a biological ancestral posterior.
Long-run construction, native custody/admission, memory/CPU/file bounds and
aggregate disk monitoring remain required. The measured
[native-plus-array disk scenario](ancestral-derived-horizon-costs-20261006.md)
remains approximately 1.44 TB before metadata growth and temporary products.
All eight evolutionary aims remain incomplete. No new MCMC, GPU prediction
or paid infrastructure is launched here.
