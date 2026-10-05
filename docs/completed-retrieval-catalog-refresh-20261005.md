# Completed retrieval snapshot and full catalog refresh

The original retrieval queue reports completion of all 1,589,879 queued
accession attempts on October 5. Individual attempts can fail or return no
exact full-length model. Queue completion does not establish complete
structural coverage of the study.

The full log was frozen under the original exclusive retrieval lock after
checking that the exact original wrapper was absent and no retrieval script
was running. All 3,538,619,125 bytes were copied with no incomplete tail;
source and snapshot hashes agree. Streaming checks cover all 3,020,697
records: 3,010,575 verified records, 9,844 errors and 278 no-exact-model
records. These are append-record counts, not distinct current accessions,
proteins, taxa or scientific coverage.

[Snapshot validation](../metadata/completed_afdb_inventory_freeze_20261005_v2.json)
and [original snapshot execution closure](../metadata/completed_afdb_inventory_freeze_transport_20261005_v2.json)
retain the native/wrapper identities, actual original tool terminal,
complete invocation journal, source hashes and declared resources.

The first snapshot adapter failed before creating its output root because
it required retained systemd metadata for the completed retrieval service.
That transient service was already missing. Its default inactive/success
fields are not proof of an exit status. The failed V1 code, logs, tool result
and [failure audit](../metadata/completed_afdb_inventory_freeze_failure_audit_20261005_v1.json)
remain immutable. The new V2 snapshot uses process absence, the completed
queue receipt/log, matching manager invocation and exclusive locking.
The original retrieval native exit status is unavailable; the snapshot
does not invent one or restart retrieval.

## Full atlas update now running

The new [full catalog plan](../metadata/completed_afdb_catalog_plan_20261005_v1.json)
retains **all 526 entries: 501 fungi and 25 outgroups**, and every one of the
**5,815,847 representative proteins**. The qualified catalog implementation
is unchanged. It uses the latest record for each accession, excludes latest
failed/no-match accessions, selects GDM/AlphaFold Monomer v2.0 pipeline
models by exact sequence hash and length, and ranks candidates by mean
C-alpha pLDDT, version and model identifier. All selected coordinate files
are rehashed against their verification records.

The new adapter binds the immutable raw producer receipt and output hashes.
An independent reader is queued behind the exact original producer. Its
separate FASTA parser, sorted ranking, every protein link, selected-model
provenance and per-taxon counts are preserved exactly from the qualified
reader. Its replacement execution gate checks the actual original
invocation, full wrapper header/terminal payloads, native exit and receipt
hash. It does not rely on systemd success defaults. The reader does not
repeat coordinate byte hashing or add a CIF-content/confidence audit.
[Source equivalence review](../metadata/completed_afdb_catalog_source_equivalence_20261005_v1.json).

The producer has original tool session **20607** and the queued reader has
session **32067**. Their immutable launch records are
[producer](../metadata/completed_afdb_catalog_launch_20261005_v1.json) and
[reader](../metadata/completed_afdb_catalog_readback_launch_20261005_v1.json).
At 15:29 UTC, 50,000 of 2,935,733 selected coordinate files have logged hash
verification. The reader is still waiting; completed coverage is unverified.
[Original process/cgroup checkpoint](../metadata/completed_afdb_catalog_checkpoint_20261005_goal_1531.json).

Each active stage has two CPU equivalents, a 96 GiB/no-swap cgroup,
88 GiB address-space cap and one BLAS thread. The producer has a 16 GiB
output allowance and a 0.5–36 hour uncalibrated planning range. The previous
producer took 2,229.93 seconds for 723,178,056,701 coordinate bytes; this
larger snapshot and disk contention make that an uncertain comparison.
The reader has a 1 GiB output allowance and a separate 0.1–12 hour
uncalibrated stage range. Its historical elapsed time includes waiting and
is not a standalone reader benchmark. Stage CPU and wall limits are caps,
not ETAs. There are no new downloads, predictions, GPUs or charges.

## Reproduction and outstanding work

Raw inventory and catalog tables remain outside Git:

- `results/structures/afdb-completed-inventory-20261005-v2/`
- `results/structures/whole-proteome-afdb-catalog-20261005-v1/`

The snapshot SHA-256 is
`aa3b3e770e1ae46b6ea3630d08722f8e0ad3a7b845effe71e94ce3a9381c4eb1`.
The plan pins representative FASTAs through their existing inventory,
the complete sampling manifest, prediction source, original scripts and
snapshot receipts. Raw model JSONL retains accession, version, file path,
sequence and coordinate hashes, software/provider and confidence metrics.
Public release of the large inputs remains outstanding; cloning the
repository alone does not supply them. The existing acquisition workflow
is `scripts/retrieve_matched_models.py` with source-specific provenance
preserved in its raw log.

Reproduction requires those immutable inputs and fresh output paths in new
pinned plans. The original producer command is:

```bash
python scripts/run_completed_afdb_catalog_v1.py \
  --plan metadata/completed_afdb_catalog_plan_20261005_v1.json \
  --receipt metadata/completed_afdb_catalog_20261005_v1.json
```

For a new run, first replace the output/receipt paths and rebuild input
pins; do not rerun this command over an existing namespace. The independent
reader requires a launch record matching the new original producer and
its completed execution receipt. Observe the current original jobs without
restarting them using:

```bash
python scripts/record_completed_afdb_catalog_checkpoint_v1.py \
  --output metadata/completed_afdb_catalog_checkpoint_NEW.json
```

After both actual original waits and the full reader close, compare the
new catalog to the September 28 snapshot across every taxon and sequence.
Then refresh confidence/PAE accounting, explicitly integrate ESMFold,
and build domain/structural-family analyses with orthology and sequence
controls. Snapshot availability is not confidence-qualified atlas
completion, homology, an accepted phylogeny or an evolutionary result.
All eight scientific aims remain required and incomplete.
