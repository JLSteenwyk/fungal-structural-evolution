# Full historical BAli-Phy scalar and rate assessment

Read-only assessment and separate serialized readback now cover all 1,620
historical short-sampler roles, 405 model groups, 135 effective inputs and
324 original aliases. The assessment found **14 rate frames with numeric
discrepancies**, including **two that pass the strict mean-one check**. Scalar
JSON also contains **33 altered alpha values** and **19 malformed records**.
These findings extend the earlier twelve-frame normalization census; they do
not repair historical output or establish adequate ancestral posteriors.

| Full checked scope | Result |
| --- | --- |
| Original scalar rows assessed | 34,020 |
| Valid scalar JSON rows compared with mapped TSV | 34,001 |
| Scalar numbers compared | 1,462,043 |
| Observed numeric scalar discrepancies | 33, all in ASRV.Gamma:alpha |
| Malformed scalar JSON rows | 19, retaining invalid bare `inf` tokens |
| Numeric comparisons unavailable in malformed rows | 817 |
| Original property frames | 4,860 |
| Category rate values compared | 19,440 |
| Discrepant category values | 17 in 14 frames |
| Corresponding repeated amino-acid rate cells | 340 |
| Discrepant frames passing strict normalization | 2 |
| Historical rate cells remaining unqualified | All 388,800 |

Each category rate repeats across twenty amino-acid states; 340 repeated cells
are not 340 independent observations. The two additional frames are iteration
zero of broad-prior roles for OG0001200 and OG0000107. Their first rates are
logged as approximately 2.022e-29 and 1.007e-21; the pure native references at
their TSV alpha values are approximately 2.022e-290 and 1.007e-210. Both changes
are enormous relative errors but negligible in the sum of four rates. That is
why normalization alone missed them. The complete discrepancy ledger retains
original values and source hashes, without replacing either number.

## Source and reference checks

The original sampler's full output/two-journal closure is inherited. All
consumed plans, role identities, native receipts, column mappings, TSV/JSON
logs and property files are freshly hashed against that closure. Scalar JSON
records are flattened by the producer and resolved directly through their
mapped hierarchical keys by a separate reader. All numeric fields are
compared at relative tolerance 2e-13 with zero absolute tolerance. Malformed
records remain explicit; the parser does not replace `inf`, infer a damaged
exponent or manufacture missing comparisons.

A single pure installed BAli-Phy calculation evaluates the original four-bin
Gamma conditional means at every saved-frame TSV alpha. It uses the installed
`gammaRatesMean` function, applies the mixture's mean-rate scaling and emits
native CJSON values plus exact encode/read checks and the original encoder's
rendering. All 19,440 native rate roundtrips pass. Native and Python parsing of
the logged finite alpha tokens agree exactly for this dataset. This is not
proof that a decimal log recovers every original sampler binary64 value.
All fourteen discrepant frames agree with the faulty original encoder's
rendering at the declared precision.

SciPy incomplete-gamma functions provide a second computational reference for
conditional bin means; the large-alpha branch uses the normal-limit formula
shown in the pinned primary implementation. All 4,860 frames agree with native
rates at relative tolerance 1e-7 and absolute tolerance 2e-12. This broader
cross-library comparison is distinct from the zero-absolute-tolerance
historical discrepancy test. It does not detect arbitrarily tiny serialization
damage, establish gamma-model adequacy or qualify the original historical rates.
The reader shares this SciPy reference function, while its scalar traversal is
separate; no third independent numerical implementation is claimed.

Primary references are the pinned
[native distribution source](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/src/builtins/Distribution.cc),
its installed ASRV/mixture/logger modules, and SciPy 1.15.3
[inverse incomplete gamma](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.special.gammaincinv.html)
and [incomplete gamma](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.special.gammainc.html)
documentation. Source copies and installed code hashes are in the assessment
bindings. Mathematical reference calculations do not substitute reconstructed
values for damaged historical output.

## Execution, failure preservation and reproduction

The first producer completed all source scanning and its pure native batch,
then failed because its parser treated the native timing footer as JSON.
Original wait 69576 exited one. The failed script, inputs, full native output,
wrapper journal and receipt are preserved in their V1 namespace. The native
child itself exited zero; no MCMC was run. Fresh V2 recognizes exactly 4,860
JSON records followed by the complete six-line footer, without changing any
numeric tolerance or historical source.

V2 original wait 88091 and separate readback wait 75257 both exited zero.
The exact original wrapper PID/create/command, invocation-specific manager
start/end and entire terminal receipt/hash payload were verified for each.
Full closure binds 8,170 source/artifact paths and both original driver journals.
Producer CPU/wall/reported child peak RSS were 34.36 seconds/34.49 seconds/
346,537,984 bytes; readback values were 21.08 seconds/21.15 seconds/
216,764,416 bytes. Reported child RSS is distinct from a cgroup peak.

Both drivers used two CPUs, 16 GiB, no swap, one BLAS thread, 12 GiB address
space, 3,600 CPU seconds, 7,200 wall seconds and 256 MiB per file. The pure
native child used 8 GiB address space, 600 CPU seconds, 900 wall seconds and
64 MiB per file. Minimum free RAM/disk were 16/128 GiB; planning allowed
8 GiB output. No GPU, new charge, MCMC, native retry, installed patch or
historical-output edit occurred.

Large ledgers and original pure native files remain outside Git at
`results/ancestral/full-historical-number-assessment-20261004-v2/`. Repository
receipts include:

- `metadata/full_historical_baliphy_number_assessment_20261004_v2.json`
- `metadata/full_historical_baliphy_number_readback_20261004_v1.json`
- `metadata/full_historical_baliphy_number_assessment_completed_20261004_v2.json`
- Both execution/transport records and the preserved V1 failure proof.

Use fresh destinations with the declared resource limits:

```bash
python scripts/assess_full_historical_baliphy_numbers_v2.py \
  --output NEW_ASSESSMENT_ROOT --receipt NEW_ASSESSMENT_RECEIPT.json
python scripts/readback_full_historical_baliphy_numbers.py \
  --producer NEW_ASSESSMENT_RECEIPT.json --transport ORIGINAL_VERIFIED_TRANSPORT.json \
  --output-root NEW_ASSESSMENT_ROOT --receipt NEW_READBACK_RECEIPT.json
```

The transport must be created from the actual original bounded tool wait and
invocation journal; a manually invented success receipt is insufficient.

## Remaining correction and scientific work

The V5 logger fixes ancestral-property encoding, while scalar parameter JSON
still uses the faulty path. Before future inference relies on it, prepare and
qualify a separate reversible source version that also emits scalar parameters
correctly, preserves existing TSV/model/distribution/transition behavior, and
records nonfinite limit events explicitly rather than disguising them as
ordinary finite values. Keep all current programs and attempts immutable.
The installed formatter remains unfixed. This finite set of observed
discrepancies is not a complete census of unobservable historical float loss.

Full posterior/mixing/likelihood/root/model uncertainty and the historical
crashes remain unresolved. The ancestral structural predictions, mechanistic
case studies, complete structural atlas, accepted phylogenetic/reconciliation/
dating framework and all eight biological aims remain incomplete. The original
weighted covariance and timing jobs continue; GPU prediction stays paused.
