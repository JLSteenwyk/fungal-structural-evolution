# Full atlas PAE inventory, validation and retrieval — October 5, 2026

## Update: full independent reconstruction closed at 20:31 UTC

The complete independent decoder/statistic reader now passes all **55,959
matrices, 112 jobs and 18,290,462,939 directional entries**. All numeric digests,
identities, bounds, extrema and sorted threshold counts match. The largest
absolute mean difference is 3.552713678800501e-15, at most 0.005997 of the
predeclared roundoff bound. Both sources remain separate; no symmetrization,
changed source bounds or evolutionary tolerances occur. The original reader
API session 81310 has an actual zero terminal and complete native/full-journal
closure, with all 61,660 declared bindings reverified. This closes independent
available-matrix reconstruction; the original producer API terminal remains
unavailable as recorded. Full missing-matrix acquisition, native context
qualification and accuracy calibration remain outstanding.

- [Complete independent reconstruction](../metadata/full_atlas_pae_statistics_readback_20261005_v1.json).
- [Original reader execution closure](../metadata/full_atlas_pae_statistics_readback_transport_20261005_v1.json).

At 20:39 UTC the original corrected downloader has **4,900 verified retrievals
and zero failures**; its complete 2,905,096 queue and full independent cache/body
readback remain unfinished. Original coordinate auditing reaches 640,000
models/247,159,346 positions; domain extraction reaches 689,242 intervals; full
timing reaches 1,248 checkpoints with 16 workers. Full biological fits remain
zero. The older running-reader checkpoints below are preserved as history.

This stage retains all 2,961,055 original source models across the fixed
501-fungus/25-outgroup representative universe. Predicted aligned error (PAE)
matrices support later assessment of relative residue/domain placement
uncertainty. Retrieving PAE does **not** add a protein structure or establish
prediction accuracy, physical domain boundaries, homology or evolutionary events.

## Complete frozen inventory and independent availability reconstruction

The original cache was frozen under its exclusive retrieval mutex. All 30,694
original receipt/matrix pairs were byte-bound; 57 lie outside the current model
selection and remain explicit. Every original model, sequence, length, version
and advertised URL is retained. The independent classifier imports no producer
classification code and reconstructs every model disposition and the entire
ordered missing queue. Both original native/API executions and complete journals
are closed, with all declared source hashes verified.

| Source | Models in inventory | Available matrix candidates | Missing advertised matrices |
| --- | ---: | ---: | ---: |
| AFDB | 2,935,733 | 30,637 | 2,905,096 |
| ESMFold | 25,322 | 25,322 | 0 |

There are no cache provenance conflicts or missing advertised URLs in this
snapshot. Existing candidates contain 18,290,462,939 directional matrix entries;
the missing AFDB queue represents 642,222,712,921 entries. Counts include diagonal
and both directional entries, which are not independent biological observations.
The unlinked annotated ESMFold gene product and all predictor configurations are
retained. Availability is separate from acceptance for scientific analyses.

- [Frozen full inventory](../metadata/full_atlas_pae_inventory_20261005_v1.json)
  and [original execution closure](../metadata/full_atlas_pae_inventory_transport_20261005_v1.json).
- [Every-model/cache/queue independent reconstruction](../metadata/full_atlas_pae_availability_readback_20261005_v1.json)
  and [original execution closure](../metadata/full_atlas_pae_availability_readback_transport_20261005_v1.json).

## Validation of every available original matrix

Four CPU workers validated all **55,959** original available matrices: 30,637
exact AFDB matrices and 25,322 ESMFold matrices. The stage replays the complete
source inventory and missing queue before decoding available bodies. It checks
source hashes/receipts, sequence/version/length provenance, complete square
dimensions, finite nonnegative values and unchanged source-specific export
bounds. AFDB retains the existing +0.51 export-rounding allowance; ESMFold retains
its existing +1e-4 allowance. Directional values and original numeric dtypes are
hashed without symmetrization, rounding or diagonal removal. Extrema, means and
inclusive 5/10/20 Å counts are descriptive summaries, not confidence admission.

All 112 jobs finished; the original native process exited zero at 19:16 UTC
after 682.99 wall seconds and 2,341.92 CPU seconds. The actual recorded native
peak RSS was 218,460,160 bytes. Source matrices occupy 29,950,413,333 bytes.

During a later executor failure, the original API wait handle became unavailable.
The immutable launch, configuration, native process receipt and complete original
wrapper journal now independently match, including the native zero exit and
original manager start/completion. **The original API terminal is unavailable**;
it is not replaced with a fabricated zero exit or a replacement command. This
native evidence is distinct from the original API-wait transport contract.

- [Full matrix validation](../metadata/full_atlas_pae_matrix_validation_20261005_v1.json).
- [Recovered original native evidence](../metadata/full_atlas_pae_matrix_validation_recovered_native_evidence_20261005_v1.json).
- [Full independent statistic reader](../scripts/readback_full_atlas_pae_matrix_statistics_v1.py)
  and [literal corruption controls](../metadata/full_atlas_pae_statistics_fixture_20261005_v1.json).

A separate full independent reader is running with four CPUs/32 GiB/no swap,
original API session 81310. It decodes original matrices without importing the
producer decoder, reconstructs every numeric digest and summary, uses sorted
order statistics for counts and extended precision for means. The declared
positive-sum float64 roundoff bound is
`eps64 * (8 * ceil(log2(max(2, N))) + 16) * max(1, abs(mean))`.
PAE value bounds and all pre-existing evolutionary fitting tolerances are
unchanged. Literal controls reject altered counts, means, numeric hashes,
symmetrization and dtypes. At 20:22 UTC, 14/112 jobs and 7,000 AFDB matrices are
reconstructed; **full reader completion and its original execution closure remain
pending**. Shared NumPy, compression and NPZ libraries remain explicit dependencies.

## Original retrieval failure and corrected complete-queue run

The first complete-queue downloader passed a relative cache directory to an
unchanged retriever that computes paths relative to an absolute repository root.
It downloaded/wrote compressed matrices, then raised `ValueError` before writing
receipts. All 8,400 recoverable error rows show this path failure; the exact
original live unit was stopped after its PID/creation/command/invocation matched.
Its immutable sources and outputs were preserved.

The original gzip lacks its end-of-stream trailer after stopping. A separate
bounded forensic review preserves that original gzip unchanged and recovers only
complete JSON lines into a new artifact. It matches all recovered rows against
the original bounded queue prefix and inventories every unreceipted file:
**8,464 files, 548,209,214 bytes, zero accepted matrices**. No HTTP provenance,
matrix acceptance or original native zero exit is inferred. Original unreceipted
bytes are not adopted or overwritten.

- [Exact original stop evidence](../metadata/full_atlas_missing_pae_stop_review_20261005_v1.json).
- [Stopped-output prefix and byte inventory review](../metadata/full_atlas_stopped_pae_review_20261005_v1.json).
- [Offline path regression controls](../metadata/full_atlas_pae_cache_path_fixture_20261005_v2.json).

The corrected V2 changes only cache-directory resolution before calling the
unchanged retriever. Offline synthetic responses reproduce the original failure
and check successful receipt writing, exact provenance, cached revalidation and
unchanged matrix rules. An initial fixture launch omitted the working directory
and failed before the native wrapper started; that failure is retained. A
distinct launch with explicit working directory passes with its actual original
tool/native zero exits and complete wrapper journal.

V2 runs the **entire original 2,905,096-model queue** with two HTTP workers and a
fixed 64-future window in separate new output/cache roots. This is an explicit
corrected full run, not an automatic retry or a corpus pilot. All successes/errors
retain original model records and receipt/matrix hashes. The first verified
20:22 UTC checkpoint has **1,100 retrieved/verified, zero failures**. These are
retrieval checks; full independent cache/queue/body readback and a fresh union
remain necessary. There are no new structures, GPU jobs, paid services or charges.

## Resources, provenance and remaining work

The missing queue's dimensional totals and existing cache imply approximately
242 GB of compressed cache and 1.70 TB of JSON transfer, with substantial
compression/length/network uncertainty. The run uses two CPUs/16 GiB/no swap,
two HTTP workers, a 2 TiB cache ceiling, a 2 TiB initial free-space requirement,
and a 1 TiB emergency reserve. More than 409 million free inodes were checked.
The 4–120 day planning range is uncalibrated; the 120-day wall/CPU cap is not a
finish ETA. The independent reader has a 48-hour wall cap and a 0.2–24 hour
uncalibrated range informed by complete original validation timings.

- [Corrected full retrieval plan](../metadata/full_atlas_missing_pae_plan_20261005_v2.json),
  [resources](../metadata/full_atlas_missing_pae_resources_20261005_v2.json),
  [original launch](../metadata/full_atlas_missing_pae_launch_20261005_v2.json)
  and [checkpoint](../metadata/full_atlas_missing_pae_checkpoint_20261005_goal_recovered_2025.json).
- [Full statistic reader plan](../metadata/full_atlas_pae_statistics_readback_plan_20261005_v1.json),
  [resources](../metadata/full_atlas_pae_statistics_readback_resources_20261005_v1.json),
  [original launch](../metadata/full_atlas_pae_statistics_readback_launch_20261005_v1.json)
  and [checkpoint](../metadata/full_atlas_pae_statistics_readback_checkpoint_20261005_goal_recovered_2025.json).

Plans and receipts bind immutable scripts, inputs and artifacts with SHA-256.
Large cache files, complete inventories, jobs, reports and interrupted byte
manifests remain outside Git under their recorded `data/` and `results/` roots.
Execution configuration records the exact reproducible bounded command; scripts
reject existing output roots. Re-execution requires separate versioned plans and
roots, rather than overwriting retained evidence.

## Queued full missing-queue independent readback

The corrected complete-queue downloader remains the only process that can add
new PAE cache files. Its subsequent independent reader is already queued as
`fungal-full-atlas-missing-pae-readback-20261006-v1.service`, with an immutable
[wait plan](../metadata/full_atlas_missing_pae_readback_wait_plan_20261006_v1.json)
and [launch record](../metadata/full_atlas_missing_pae_readback_launch_20261006_v1.json).
The waiter pins the exact live downloader launch, producer plan, resource
contract, reader and controller sources. It checks the downloader's original
PID/creation time/command while live, then requires original systemd completion
resource evidence and a successful terminal state before the readback begins.

The reader independently replays all 2,905,096 ordered queue rows. For every
successful retrieval it compares the queued source identity with the retained
receipt, verifies receipt/compressed/decompressed SHA-256 values and byte counts,
then checks the JSON container, square dimensions, finite/nonnegative entries
and the established AFDB export bound. Retrieval failures remain explicit and
their error-type census must match the producer. It does not call the downloader
or import its PAE validator. Literal malformed-container/shape/value/hash
controls and a real retrieved receipt/matrix check pass before queueing.

The waiting and reader stage is limited to **2 CPU equivalents, 16 GiB memory,
zero swap and no GPU**. It has a 12 GiB address-space cap, three-day CPU/wall
limits and a fresh compact output root. The 4--72-hour readback range is
uncalibrated. This check does not symmetrize PAE, determine confidence or
domain-context eligibility, assess predictive accuracy, add structures, or
support biological inference.

### Live completed-prefix checkpoint

Before the complete downloader finishes, a separate read-only checkpoint
validated its first 1,000 completed records. Completion order differs from
source-queue order because two HTTP workers finish asynchronously, so the reader
maps every observed prefix model back to its unique original queue record before
checking it. All 1,000 records passed source identity, retained receipt,
compressed/decompressed SHA-256 and byte-count checks. The checkpoint also
decoded **56,250,382 compressed bytes**, **419,607,757 JSON bytes** and
**169,326,412 directional PAE values**, requiring square dimensions and finite,
nonnegative values within the established AFDB export bound.

The [checkpoint](../metadata/full_atlas_missing_pae_prefix_checkpoint_20261006_v1.json)
and its [independent reader](../scripts/check_live_missing_pae_prefix_v1.py)
do not read the changing end of the live gzip stream, accept the remaining
queue, or replace the required whole-queue post-completion reader.

Whole-source coordinate auditing and full domain extraction remain running:
at 20:20 UTC, 562,000 models/214,165,994 residues are audited and 610,020 domain
intervals exported. Full timing has 1,231 producer checkpoints and 16 actual
workers; full weighted biological fits remain zero. These are partial runtime
observations, not completed independent scientific analyses. GPU prediction
remains paused. Full structural coverage, context/accuracy calibration,
source-aware comparisons, accepted phylogenetic/reconciliation/dating
frameworks, adequate ancestral uncertainty, all eight evolutionary aims and
final data/figure/methods/case-study deliverables remain unfinished.
