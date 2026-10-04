# Scalar JSON V6 qualification and future source grid

The project-local V6c logger corrects future scalar parameter serialization,
retains explicit nonfinite states, and avoids the installed native literal-null
encoder. It is software qualified and prepared across all 405 V5 programs and
1,620 future chain roles. It has **not** been integrated into the production
sampler audits or launched across the full native grid. Historical numbers,
installed libraries, existing model sources and current jobs remain unchanged.

## Representation and evaluation

The probability model, priors, reference-alignment initializer, original TSV
logger and joint ancestral logger are preserved byte-exact by reversing the
source transformation. Scalar context fields are evaluated once with the same
`runContextAction`; the resulting values receive quality metadata before native
encoding. No extra probability action, ancestral draw or logger is registered.
Parameter fields pass through the native CJSON finite-number encoder instead
of the faulty Haskell Double formatter.

The MCON header declares `projectScalarSchema` as
`native-cjson-explicit-special-values-v6`. Each row contains `iter`,
`statistics//`, `parameters//`, and `numericParameterQuality//`. The context
object adds `__project_scalar_v6_quality__`. Both quality objects carry
`numericLeafCount`, `nonfinite`, and `literalNullPaths`. Paths retain string
keys and integer array indices as separate typed segments.

Nonfinite values become reserved strings:
`__project_scalar_v6__:positive_infinity`, `negative_infinity`, or `nan` with
the same prefix. Literal nulls become `__project_scalar_v6__:literal_null`.
Their paths and kinds must match the quality records exactly. These strings
are diagnostic values, never numerical substitutes or qualified samples.
Reserved context keys and original strings with the reserved prefix are
rejected by the helper. Production consumers must explicitly support this
schema; the existing sampler audits are incompatible and remain unchanged.

Two separate Python validators traverse the records independently. They reject
unrepresented nulls/nonfinite values, duplicate or overlapping tags, malformed
typed paths, wrong placeholder kinds and incorrect numeric counts. The separate
reader also maps actual finite model records onto original TSV columns at
relative tolerance 2e-13, absolute tolerance zero. It refuses diagnostic
nonfinite rows for that finite-only comparison; it does not accept them for
posterior inference.

## Actual qualification

The successful V9 checker qualified all 405 source transformations and all
1,620 role identities statically. Three V5/V6 paired native comparisons use
the same fresh fixture seeds for broad, centered and package priors, each for
20 iterations. This is software verification, not a separate biological pilot.

All original TSV files, column maps, trees, ancestral FASTA files and joint
site-property JSONL files were byte-identical across each pair. Both validators
checked 63 corrected scalar rows; independent mapped TSV comparisons checked
2,709 numeric values. Both existing joint-frame decoders and serialized NPZ
readback checked nine frames and 36 ancestors without changing normalization
or sequence/state/coordinate assertions.

A separate pure-native Value fixture covered 18 finite values, extreme
exponents, the smallest positive subnormal and signed zero, plus six explicit
nonfinite tags and two literal-null tags across context and parameter sections.
Seventeen finite values passed exact native encoding/readback, including zero
sign verified by reciprocal sign and Python `copysign`. The installed
`Text.Read` refuses `5E-324`; the subnormal remains tested separately using
native `encodeFloat`/`decodeFloat`, native components
`[4503599627370496, -1126]`, exact Python `ldexp` reconstruction and exact JSON
readback. No native Text.Read subnormal roundtrip is claimed. Twenty-four
altered records were rejected by both independent validators.

The actual original tool wait was **89239, exit 0**, wrapper PID **1301056**,
invocation **73f9ebfde4a94fc781ca17f0d010182e**. Both complete wrapper messages
and original manager start/end were verified, including the entire terminal
receipt/hash payload. The transport proof freshly checks 4,225 bindings.
The stage used 2 CPUs, 16 GiB RAM, zero swap, one BLAS thread, AS 12 GiB,
3,600 CPU seconds, 7,200 wall seconds and a 256 MiB per-file ceiling.
Observed child CPU/wall time was 58.23/58.41 seconds; peak child RSS was
808,341,504 bytes, not a cgroup peak. The shared wrapper's inherited scope text
names an older kernel checker; its actual command and the transport proof
identify this scalar checker unambiguously.

## Preserved failures and native limitations

All eight unsuccessful software namespaces remain intact with their original
actual waits, complete wrapper payloads and manager journals. None was restarted
or overwritten; corrections used new checker versions and new roots.

| Checker | Retained failure |
| --- | --- |
| V1 | Fixture name `values` collided with an installed import. |
| V2–V4 | Native `float_next(inf)` domain error with extreme fixture construction; replacing the maximum finite literal alone did not resolve it. |
| V5 | Pure Value fixture exited SIGSEGV. |
| V6 | Correct special-value JSON was written; installed Text.Read rejected the smallest subnormal. |
| V7 | Native output completed; the checker incorrectly assumed a mantissa convention and trusted `isNegativeZero`. |
| V8 | Native compilation rejected the unavailable `ToJSON Integer` instance. |

The four-program diagnostic's actual wait **16702, exit 0** and original
journal are independently retained. GDB localized the V5 crash to
`builtin_function_ejson_null` via `OperationArgs::reg_for_code`. GDB's zero
exit is **not** native success. Ordinary finite quality encoding, extreme
finite encoding and primitive nonfinite inspection exited zero. Primitive
native nonfinite output was `[Infinity,-Infinity,NaN]`, which is not strict
JSON. Thus native CJSON alone cannot provide nonfinite-safe scalar logging.
This evidence corrects the earlier assumption that it would emit null safely.
Sanitization avoids both observed paths without patching installed software.

## Complete prepared grid and next work

The prepared manifest retains 405 model quartets, 135 effective inputs,
324 original aliases, and 540 roles per prior. All 1,620 new seeds are unique,
order-independent and disjoint from a 6,518-seed census of original sampler,
earlier future, stack-followup and prior BAli-Phy software namespaces. Fixture
seed reuse exists only within declared paired software comparisons. The future
receipt binds 4,634 source/artifact paths and excludes fixture constants/oracles.
No production horizon or execution command was prepared.

Next: integrate a versioned V6-aware full sampler audit, qualify the full native
source/runtime grid under a documented resource plan, and resolve the largest
family failures before selecting posterior horizons. Independently diagnose
the retained covariance timing failure; its original 10 partial cohort
receipts and all three failed controllers are preserved. Neither original
numerical tolerances nor prior outputs may be relaxed or overwritten.
Weighted numerical qualification continues independently. Real full-data fits,
calibration, complete atlas/framework and all eight biological aims remain open.

## Reproducible artifacts

- Helper: `config/native-format-probes/project-scalar-json-v6c.hs`.
- Transformation/primary validator: `scripts/baliphy_scalar_json_logger_v6c.py`.
- Independent reader: `scripts/read_baliphy_scalar_json_v6b.py`.
- Checker: `scripts/check_baliphy_scalar_json_logger_v6_v9.py`.
- Validation: `metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json`.
- Execution/transport: matching `...software_execution...v9.json` and `...software_transport...v9.json`.
- Diagnostic: `metadata/baliphy_scalar_json_v6_diagnostic_20261004_v1.json` and matching transport.
- Prepared grid: `metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json`.
- Timing failure: `metadata/full_retained_shared_entity_timing_failure_20261004_v2.json`.

Large generated/native artifacts stay outside Git under the bound
`data/software_audits/` and `results/ancestral/` paths. Receipts preserve their
exact locations, configurations, process identities and hashes.
