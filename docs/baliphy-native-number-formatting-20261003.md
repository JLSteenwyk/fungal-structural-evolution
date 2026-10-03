# Native scientific-number formatting failure

All 1,620 historical 20-iteration sampler runs completed native execution and
output/readback closure with zero allocation failures. This is a computational
qualification result, not an adequate posterior. Their resource observer also
closed. Subsequent independent replay exposed a separate **numeric serialization
bug** in the installed BAli-Phy 4.3 Haskell JSON path.

The strict replay checks failed because logged category rates did not average
to one. A complete census of all 4,860 saved property frames found twelve such
failures across OG0000294, OG0000230 and OG0000972. The largest mean-rate
discrepancy was about 0.200. This threshold census is not a complete assessment
of historical floating-point corruption: some altered small values can remain
below the mean-check tolerance.

Primary code at commit `80b0402eed0157f31ecb57e0efc34c03ed83050c` was retrieved
and checked against Git blob hashes. Its
[Prelude formatter](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/src/builtins/Prelude.cc#L353)
removes trailing zeros from any formatted number containing a decimal point,
including zeros belonging to a scientific exponent. The installed JSON double
encoder calls that formatter. A separate pure computation with the actual
installed binary reproduced the problem: `2.34e-10` became `0.234`,
`2.34e-20` became `0.0234`, and `2.34e+20` became `234`. Twelve fixed constants
included five incorrect encodings and unchanged controls. No MCMC was run.

An independent Gamma-category calculation from matching scalar-log alpha values
also identifies this pattern in the saved records. For example, the first
problematic frame logs `0.406972398064201` where the category mean is approximately
`4.0697239806419846e-10`; its remaining three categories agree with the expected
rule. This comparison supports the serialization diagnosis. It is not adoption
of reconstructed rates as an exact replacement for historical native values.
Encoded numbers alone do not uniquely identify their original exponents.

The six original replay/joint producer, reader and closure services stopped.
Their invocation journals and original sources/output directories remain
preserved. **The joint sampler started zero native roles.** Dependent services
correctly did not proceed without successful replay closure; some reported
missing completion-resource evidence after the failed prerequisite. None was
restarted and the rate-mean assertion was not weakened.
[Complete failure evidence](../metadata/baliphy_native_exponent_formatting_failure_20261003_v1.json)
and [full frame census](../metadata/short_sampler_native_rate_mean_census_20261003_v2.json).

## Tested alternative and remaining work

The installed native CJSON path, `cjsonToText . toCJSON`, preserves scientific
exponents. A separate probe passed exact **native encode/read roundtrips for all
twelve fixed values**, with no installation edits. A first CJSON checker is
preserved as failed: Python and the native compiler parsed one decimal literal
to adjacent binary64 values, differing by one ULP. The corrected checker keeps
that cross-language difference explicit and uses native encode/read equality as
its exact oracle. It also rejects the large exponent changes through the
cross-language numeric comparison. This does not loosen the sampler's rate
normalization requirement.
[Successful native CJSON probe](../metadata/baliphy_native_double_format_cjson_probe_20261003_v2.json)
and [preserved literal-oracle failure](../metadata/baliphy_native_double_format_cjson_v1_development_failure_20261003.json).

Before new joint sampling, implement a separate reversible logger version that
uses the correct numeric encoding. Qualify it across all 405 prepared programs
and 1,620 roles, with saved native frames across all priors and complete independent
sequence/category/coordinate checks. Update future plans and resource observation
to the new versions rather than changing failed sources or resetting old counters.
Historical replay must preserve every role and affected property record, document
the broader formatting assessment, and distinguish unavailable/flagged historical
rates from usable sequence/category observations. Do not infer exact original
numbers merely by renormalizing the damaged JSON.

The pure probes used two CPUs/8 GiB/no cgroup swap, with a 20-second CPU bound,
30-second wall bound and 16 MiB per-file limit for the native invocation. The
original formatter and corrected CJSON checks used about 0.90 and 1.03 child CPU
seconds. The original binary and installed libraries remain unchanged. No GPU
prediction, new infrastructure charge or existing native-chain retry occurred.
Longer posterior adequacy and all eight biological aims remain incomplete.
