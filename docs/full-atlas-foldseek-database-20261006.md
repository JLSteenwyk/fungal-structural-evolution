# Protein-level structural-search database across the expanded atlas

The full current atlas contains **2,961,055 source models**: 2,935,733 selected
AlphaFold database models and all 25,322 completed ESMFold models. The earlier
completed protein database contains 1,290,278 AlphaFold models from the September
22 catalog. Its existing clusters remain tied to that older input set. A new
database is now being built across the entire current atlas, preserving both
predictors and the one ESMFold alternative without a representative-protein link.
This adds search infrastructure for existing structures, not new predictions.

| Input quantity | AFDB | ESMFold | Total |
| --- | ---: | ---: | ---: |
| Source models | 2,935,733 | 25,322 | 2,961,055 |
| Inventory residue positions | 1,152,706,482 | 12,221,520 | 1,164,928,002 |

These are source models and positions, not independent proteins, accepted
homologous families or biological observations. The existing availability union
retains all 5,815,847 representatives across 501 fungal entries and 25 outgroups,
with missing models and predictor overlaps explicit. No predictor is chosen
by confidence and no species-specific source link is discarded for indexing.
Entry counts still require taxonomic/assembly review before being treated as
verified unique species.

## Conversion and complete source-coordinate checks

The new builder requires the closed complete AFDB/ESMFold availability union,
the exact full model-metadata digest, the qualified installed Foldseek binary
and completed mixed-source software controls. It enrolls every input, checking
source labels, lengths, digests, CIF paths and unique filename aliases. A missing
input or colliding alias fails before native conversion; it is not an exclusion.

Native `createdb` uses the same installed executable and scientific options as
the older qualified database: **eight CPU threads**, GPU zero, pLDDT/B-factor
threshold **70 for seeding**, and coordinate store mode **1 (float32)**.
`CUDA_VISIBLE_DEVICES` is empty. Seeding masks do not establish residue-confidence
eligibility or calibrated prediction accuracy. Amino-acid sequences, native 3Di
states, coordinates and original source identity remain traceable.

After native conversion, the existing qualified database checker verifies every
lookup identity, amino-acid sequence/hash, 3Di alphabet/length, coordinate shape,
finiteness and record terminator. A separate stage then freshly hashes every
original CIF, reconstructs its sequence and atom/residue/chain/alternate-location/
occupancy/confidence properties with the already qualified source-reader helper,
and compares **every native C-alpha coordinate exactly with the original value
after float32 conversion**. The full original helper implementation is unchanged.
Neither PDB rounding nor a relaxed coordinate tolerance is introduced.

The coordinate stage has 1,000 models per job and retains every source and model
identity, native offset/record size, job digest and completed shard proof. Source
and output hashes close only after the entire source-model census agrees and
the native coordinate-file digest remains unchanged. Per-source zero counters
are explicit. Partial output and the original native execution are retained if
a stage fails; no failed source model is silently dropped.

Full independent database readback and original execution closure remain
subsequent requirements. The producer's source-coordinate check is independent
of Foldseek's native CIF conversion, but uses the project's shared CIF lexical
parser. Native **3Di values are not independently reconstructed** here. Complete
3Di accuracy/uncertainty benchmarking, source effects, confidence/PAE calibration,
structural clustering, edge validation and biological homology assessment remain
separate analyses. Structural similarity groups cannot be called orthogroups.

## Software qualification and preserved first failure

The mixed-source compatibility check converts two actual AFDB and two actual
ESMFold CIFs. All four amino-acid/3Di/native-coordinate records pass, including
exact comparisons of every C-alpha coordinate. Eleven controls reject duplicate
source aliases, unknown predictors, boolean lengths, malformed digests,
newline-injected paths, amino-acid source disagreement, altered finite coordinates,
nonfinite coordinates, truncated records, wrong terminators and changed source
CIF digests. Positive native outputs are retained outside Git with checksums.
These are software controls within the full workflow, not a taxon/family pilot.

The first V1 checker failed before native conversion because it imported a
nonexistent helper name. Its original tool return completed immediately with
actual exit one and no session ID. V1 source, logs, execution and complete
original wrapper journal remain unchanged. Fresh V2 imports the existing
`independent_arrays` helper; no scientific setting, source parser or numerical
rule changes. Original V2 check session **19668** returns actual zero and has
complete source/wrapper/manager transport closure. Its transport session **97251**
also returns actual zero. The fixture does not qualify full conversion or homology.

Evidence:

- [Preserved V1 failure](../metadata/full_atlas_foldseek_fixture_failure_20261006_v1.json).
- [Mixed-source native controls](../metadata/full_atlas_foldseek_fixture_20261006_v2.json)
  and [original execution closure](../metadata/full_atlas_foldseek_fixture_transport_20261006_v2.json).
- [Unchanged source reader](../scripts/readback_full_atlas_coordinate_profiles_v1.py)
  and [full source/database builder](../scripts/build_full_atlas_foldseek_database_v2.py).

## Resource estimate and actual full launch

The previous complete native conversion processed 1,290,278 models and
514,418,704 residues with four CPU threads in **13,455.124 seconds
(3 h 44 min 15.124 s)**. Its complete native database occupies **7,422,261,481
bytes**. Model- and residue-scaled eight-thread scenarios give **4.29 and 4.23
hours** for current native encoding, before the new full original-CIF geometry
readback. These are uncalibrated extrapolations, not a finish ETA.

The current uncompressed coordinate payload alone is **13,982,097,079 bytes**.
Amino-acid, 3Di and coordinate payloads total at least **16,323,797,303 bytes**,
before headers, indices, path/source maps and temporary products. The per-file
cap is therefore **64 GiB**, rather than inheriting a smaller limit that would
truncate the coordinate database. Prior byte counts are filesystem measurements;
the resource estimate does not claim a fresh hash audit of every older database.

The full stage has **eight CPU equivalents, 96 GiB memory, no swap, 88 GiB
address space per process and one BLAS thread**. Planning allows **250 GiB output
and 8–72 hours**; the actual wall cap is 72 hours and per-process CPU cap is
seven days. Launch requires 96 GiB available RAM and 2 TiB free disk. The estimate
observed **642.8 GiB available RAM and 9,756.7 GiB free disk**. These allowances
are not measured long-run upper bounds. No new charge or infrastructure is used.

Original full tool session **33972** launched
`fungal-full-atlas-foldseek-database-20261006-v2.service`, invocation
`12cb88ce14b54334ac3f83b6a11848b0`. Original wrapper PID **765873**, creation
time **1791261371.56**, exact command and actual eight-CPU/96-GiB/no-swap kernel
limits are recorded. At 04:37 UTC the original process was live and had enrolled
900,000 input records; native conversion and full source-coordinate completion
were not yet claimed. Poll the original session or fingerprint its original
process/invocation; an observation failure must not trigger a duplicate run.

At **04:42 UTC**, every one of the **2,961,055 models** had passed enrollment.
Source counts and residue totals match the frozen full atlas; all filename
aliases are unique and all inputs exist. Native encoding is confirmed live as
PID **768169**, creation time **1791261606.94**, using eight threads and GPU
zero with an empty CUDA-device environment. The complete source-file census is
**1,105,396,373,092 bytes**. This is a filesystem-size census; fresh hashing of
every CIF occurs in the subsequent source-coordinate stage. The
[native-phase checkpoint](../metadata/full_atlas_foldseek_checkpoint_20261006_goal_native_v2.json)
records the original live process and exact full enrollment, without claiming
native or geometry-check completion.

[Resource calculation](../metadata/full_atlas_foldseek_resource_cost_20261006_v2.json),
[bounded resources](../metadata/full_atlas_foldseek_resources_20261006_v2.json),
[immutable full plan](../metadata/full_atlas_foldseek_plan_20261006_v2.json),
[original launch](../metadata/full_atlas_foldseek_launch_20261006_v2.json)
and [launch checkpoint](../metadata/full_atlas_foldseek_checkpoint_20261006_goal_launch_v1.json)
retain the actual scope and provenance.

The [observed environment](../metadata/full_atlas_foldseek_environment_20261006_v2.json)
records Python, NumPy, Biopython, psutil, platform, byte order and the exact
installed executable digest. These are observed versions, not an archived copy
of every installed package; existing project environments remain the
reproduction starting point.

```bash
python scripts/build_full_atlas_foldseek_database_v2.py \
  --plan metadata/full_atlas_foldseek_plan_20261006_v2.json \
  --receipt metadata/full_atlas_foldseek_database_20261006_v2.json
```

Reproduction requires all pinned data/qualification evidence, fresh output paths,
the recorded bounded wrapper and kernel/BLAS limits. Existing outputs are refused.
The full model metadata, input path list, native database, geometry jobs/proofs
and all coordinates stay outside Git history with their source/output digests.
The earlier interrupted database and recovered complete database remain intact.

At 04:38 UTC the original domain reader had checked **228,536 intervals**, the
atlas audit had processed **2,706,000 models**, and the original PAE downloader
had verified **55,200 new matrices**. Original weighted timing retained
**1,582/4,340 producer cohorts**, with full biological weighted fits still zero.
These are partial process-backed observations; full independent atlas/domain/
timing closure remains pending. GPU prediction remains paused. Full missing
prediction coverage, accepted biological frameworks, calibrated branch effects,
adequate ancestral uncertainty and all eight evolutionary aims remain unfinished.
# Full mixed-source Foldseek database

After the producer completes, the required separate reader is
`scripts/readback_full_atlas_foldseek_database_v1.py`. It opens the native
lookup, amino-acid, 3Di and coordinate index files with an independent parser
and checks every source-model alias, amino-acid hash, 3Di length/alphabet and
float32 coordinate record. It runs only against a completed producer receipt
and fresh output paths. This native-record readback is separate from the
producer's original-CIF geometry validation; neither establishes structural
homology, confidence calibration or evolutionary change.

## Queued independent native-record readback

The reader is queued as
`fungal-full-atlas-foldseek-readback-20261006-v1.service` through the immutable
[wait plan](../metadata/full_atlas_foldseek_readback_wait_plan_20261006_v1.json)
and [launch record](../metadata/full_atlas_foldseek_readback_launch_20261006_v1.json).
At queue creation, it bound the live original producer launch record, its
immutable full-atlas plan, the reader/controller sources, and the declared
resource contract by SHA-256. The wait process verifies the original producer's
PID/creation time/command while it remains live. Once it terminates, the V2
dependency runner also requires that same systemd invocation to have an original
completion-resource journal record and a successful terminal state before
running the reader. A failed, replaced, or altered producer therefore cannot
trigger this readback.

The queued unit is capped at **2 CPU equivalents, 16 GiB memory and zero swap**;
it has no GPU access. The reader itself is sequential and read-only with respect
to source structures and the native database. It has a 12 GiB address-space cap,
three-day CPU/wall limits, a 512 MiB per-file cap, and fresh output/receipt paths.
It checks all 2,961,055 aliases and native records only after the producer has
created its completed receipt. Its planned 4--72-hour range is uncalibrated.
The handoff does not start conversion, structural searches, clustering, model
prediction, or any biological inference while the producer remains active.

## Preserved readback failure and corrected fresh reader

After the producer's successful completion, the queued v1 native-record reader
stopped before producing a receipt. Its retained stderr records a `MemoryError`:
the imported digest helper attempted to read a large bound database artifact in
one allocation under the reader's 12-GiB address-space cap. No model disposition,
native record, or output receipt from that failed run is accepted.

V2 streamed every bound artifact successfully, but then reached `MemoryError`
while constructing the declared 2,961,055-entry lookup/index grid under its
12-GiB address-space cap (observed child peak RSS: 12,794,757,120 bytes). Its
execution record is retained and no output from it is accepted. The unchanged
parser is relaunched with a 32-GiB cgroup / 28-GiB address-space envelope, still
with two CPUs and zero swap. This capacity change is data-supported and limited
to the fresh v3 reader; it does not alter the producer, input database, native
record checks or any scientific method.

The compact [v1 execution record](../metadata/full_atlas_foldseek_readback_execution_20261006_v1.json)
and [v2 execution record](../metadata/full_atlas_foldseek_readback_execution_20261006_v2.json)
retain the exact commands, cgroup limits, source hashes, exit dispositions and
hashes of the locally retained logs. They are failure evidence only, not
successful reader receipts.
