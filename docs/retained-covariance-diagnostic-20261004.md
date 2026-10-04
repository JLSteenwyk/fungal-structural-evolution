# Exact failed-cohort covariance investigation

The original retained-input timing stage remains failed after ten cohort
checkpoints. Its unchanged `backend_guard` rejected nine of 25 Gram entries
at `rtol=3e-9` and `atol=2e-8`, reporting maximum absolute difference 2.56917974.
Neither its original data nor its tolerances have been changed. No production
model fit has been launched.

A new bounded diagnostic reconstructs the exact eleventh cohort,
`00c9462e56f1fcdf34d9e1f94c47c072ea4d26961d46c37b18de57bbc9ddc189`,
with 22,881 observations. It uses the frozen full-source loader, reconstructs the
complete 1,200 candidate census and original selection of 40 representatives, and
calls the unchanged guard after primary likelihood construction. Every rejected
guard will be retained as diagnostic evidence. Separate retained, independent
latent and original seven-kernel principal-submatrix calculations will be
compared. The original producer has not been restarted.

The new diagnostic is live at 06:13 UTC, still verifying original sources, with
zero group records and no receipt. Original tool session 65481, wrapper 1333144
and child 1333148 match invocation 1d4b04137c474a42b2eea4aae6de0f20. Its
limits are 2 CPU/32 GiB/no swap/one BLAS thread, 24 GiB address space, 7,200 CPU
seconds, 14,400 wall seconds and 512 MiB per file. Full source checking is I/O
heavy; the snapshot records memory.max events but zero OOM/kill events. These
are resource budgets, not a completion estimate. A separate extended-precision
reader and its resource plan are prepared but have not been launched or queued.

## Completed narrower arithmetic check

While full-source verification runs, a separately bounded check examined the
residual–species and species–species **raw** Gram entries in all 300 original
audits for this cohort. For each of five trees, it grouped 6,782 distinct species
factor rows, retained their integer multiplicities and computed the species
Gram using `longdouble` arithmetic (63 mantissa bits). The underlying factor
has shape 22,881 by 301 and dtype float64 for every tree.

All 600 checked entries agree with their original saved values under the
unchanged acceptance criterion. Direct float64 calculations also agree with
the wider reference. Maximum absolute differences between the wider reference
and saved values range from 4.41e-9 to3.33e-8 across trees. This does not
reproduce the reported 2.56917974 difference in these entries and makes a simple
rounding explanation for their basic species contractions less likely. It does
not explain the original failure or establish that the complete guard passes.
Entity–species terms, projected Grams, primary construction and all 40 original
guard calls are the scope of the still-running diagnostic.

The completed check verifies 22 consumed source bindings, including ordered case
rows, factor inputs and separately closed producer receipt links. It does not
replay the entire original source archive or qualify biological input sources.
Wider floating-point arithmetic is a numerical reference, not an interval proof.
The first transport's generic scope text referred to inherited full-source
archives; for this narrow producer, no such archive exists. Separate transport
V2 explicitly records `consumed_source_bindings_only`; original V1 is retained.

Actual original wait 77536 exited 0. Wrapper 1342417/native 1342420 and
invocation 48570f2ff5074e9c9885cd9ea0d3a4b2 match the complete original
wrapper payload and manager start/end. The execution receipt gives exact
identities if this summary is used for reproduction. The check used 56.82 child
CPU seconds, 56.90 wall seconds and reported child peak RSS 2,773,012,480 bytes;
these are distinct from the manager's smaller reported memory value.

At 06:11 UTC the independent weighted numerical stage is live and reports
530/4,340 cohorts, with 533 partial audit segments and no recorded failure files.
Its complete numerical closure remains absent. No tolerance was relaxed, no
review was promoted and no original job was restarted. GPU prediction remains
paused; all eight biological aims and the complete project remain unfinished.

Evidence:

- [Exact diagnostic plan](../metadata/retained_timing_covariance_diagnostic_plan_20261004_v1.json).
- [Fresh original live-handle checkpoint](../metadata/retained_timing_covariance_diagnostic_checkpoint_20261004_0614.json).
- [Completed scoped raw reference](../metadata/retained_species_raw_diagnostic_20261004_v1.json).
- [Original wait and precise source scope](../metadata/retained_species_raw_diagnostic_transport_20261004_v2.json).
- [Prepared extended reference](../metadata/retained_timing_covariance_extended_reference_plan_20261004_v1.json).
- [Weighted numerical checkpoint](../metadata/full_weighted_covariance_qualification_execution_checkpoint_20261004_diagnostic_0612.json).

Next: poll original session 65481 without restarting it. Once terminal, preserve
its exact wait and invocation proof, inspect every selected group and the first
sorted rejection, then run the prepared independent raw reference on its closed
exports. An unexplained failure remains unresolved even if a new calculation
passes. Whole-grid source closure, timing, fitting and calibration still apply.
