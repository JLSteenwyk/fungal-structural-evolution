# Expanded full domain manifest and coordinate extraction — October 5, 2026

The current AFDB catalog supplies **2,454,565 distinct candidate intervals from
976,357 source models**, up from 1,575,294 intervals/627,567 models in the
September 28 manifest. This is approximately 56% growth in candidate intervals.
The complete independent SQL-union check passes. Coordinate extraction has
started across the entire new manifest; it is not yet complete or independently
verified. These are annotation-defined candidates, not accepted structural
boundaries, homologous families or evolutionary events.

The predecessor retains all 2,935,733 AFDB models and 2,994,868 protein links.
The existing scientific rules define which models have candidate domains;
models lacking eligible annotations remain in the full atlas and denominators.
All four original annotation-policy alternatives enter the union, and both
alignment and envelope boundaries are retained. Sequence/span identity removes
duplicate coordinate intervals reversibly, preserving every model/hit/boundary
association. Policies and alternative boundaries are not independent replicates.
ESMFold models remain in the separate complete predictor inventory; they are
not silently merged into this AFDB-only domain registry.

| Manifest quantity | Full count |
| --- | ---: |
| Eligible source models | 976,357 |
| Distinct model/hit pairs | 1,352,782 |
| Reversible boundary associations | 2,705,564 |
| Distinct sequence/span intervals | 2,454,565 |
| Residue occurrences across intervals | 383,162,558 |
| Maximum interval length | 1,103 |

Residue occurrences include overlapping intervals and are not unique residues.
The producer reconstructs every boundary association. The unchanged independent
reader separately constructs the full SQL union and checks every emitted
identity, model, sequence, path, inclusive bound, residue count and dimension.

The first original wrapper wait **70086 failed with exit 1**. Its execution
receipt argument lacked a `.json` suffix, so the wrapper created its log
directory at the same path and failed with `FileExistsError` when writing the
final report. The domain producer had already written its complete bound output.
The exact original configuration, process identity, launch, logs, whole failure
journal and actual tool payload remain unchanged. **The original native terminal
exit code is unavailable; it is not replaced with zero.**

A separate retained-output review, original wait **39729**, closes with actual
exit zero and binds the original failure and every saved artifact. It repeats
no native domain generation. The subsequent full independent SQL-union reader,
original wait **89364**, closes with actual exit zero. Both new original tool
payloads match the complete wrapper journals and manager start/completion
records. This qualifies the saved manifest for extraction while preserving
the first reporting failure. Future wrapper execution-receipt arguments include
`.json`, keeping the receipt distinct from its log directory.

Evidence and full sources are in the
[retained-output review](../metadata/completed_domain_manifest_reporting_review_20261005_v1.json),
[review execution closure](../metadata/completed_domain_manifest_reporting_review_transport_20261005_v1.json),
[independent full reader](../metadata/completed_domain_manifest_readback_20261005_v1.json)
and [reader execution closure](../metadata/completed_domain_manifest_readback_transport_20261005_v1.json).

Coordinate extraction uses the unchanged qualified all-atom exporter, with
four workers, 32 GiB memory, no swap and one BLAS thread. It attempts every
interval in 977 shards containing at most 1,000 source models each. Original
wait **64749**, unit `fungal-completed-domain-coordinates-20261005-v1.service`,
is retained and must be polled to actual terminal completion; a transient unit's
absence or partial state file is not completion evidence. The exact launch and
actual CPU/memory/swap limits are recorded in the
[launch descriptor](../metadata/completed_domain_coordinates_launch_20261005_v1.json).

Before launching, the
[immutable full extraction plan](../metadata/completed_domain_coordinates_plan_20261005_v1.json)
and [bounded resources](../metadata/completed_domain_coordinates_resources_20261005_v1.json)
estimate approximately 9.6 hours and 245 GB of tar output by scaling a prior
completed four-worker extraction of 427,255 models/1,078,592 intervals, which
took 15,201.92 seconds and produced 107.82 GB of tar files. The 8–96 hour
planning range and 750 GiB output allowance accommodate different model sizes,
system load, hashing and I/O. The 96-hour execution cap is not an ETA. Available
storage exceeds 10 TiB; launch requires 2 TiB free and retains a 1 TiB emergency
reserve. No GPU, new prediction or paid resource is used.

For each source model the exporter checks exact coordinate bytes, polymer
sequence/hash/length, atom/residue identity, finite coordinates, occupancy and
confidence bounds, one chain/model and complete C-alpha coverage. It exports
all atoms within each interval, retains sequence correspondence and confidence,
and verifies PDB rounding at the unchanged coordinate/confidence tolerances.
Every rejection and missing-backbone position remains explicit. A wrapper binds
the complete output only after native completion. Full independent archive,
atom and disposition readback is still required; it has not been queued before
the original extraction closes. Extraction is not predictor calibration,
confidence acceptance, PAE qualification or structural clustering.

Reproduce the full saved-output check with the recorded plans:

```bash
python scripts/review_completed_domain_manifest_reporting_v1.py \
  --plan metadata/completed_domain_manifest_reporting_review_plan_20261005_v1.json \
  --receipt metadata/completed_domain_manifest_reporting_review_20261005_v1.json
python scripts/readback_completed_domain_manifest_v1.py \
  --plan metadata/completed_domain_manifest_readback_plan_20261005_v1.json \
  --receipt metadata/completed_domain_manifest_readback_20261005_v1.json
python scripts/run_completed_domain_coordinates_v1.py \
  --plan metadata/completed_domain_coordinates_plan_20261005_v1.json \
  --receipt metadata/completed_domain_coordinates_20261005_v1.json
```

These commands require the pinned data, original evidence and fresh output
paths; existing outputs are deliberately refused. Native execution must use
the recorded bounded wrapper and cgroup settings from each configuration.
Large intervals, coordinates and tar shards remain outside Git history.

The [full predictor availability union](full-prediction-atlas-union-20261005.md)
has 2,961,055 source models linked to 3,019,669 representative proteins (51.9%).
About 2.80 million representatives still lack either predictor. GPU prediction
remains paused. Full confidence/PAE controls, source-aware structural comparisons,
accepted species/gene/reconciliation/dating frameworks, calibrated branch effects,
adequate ancestral posterior uncertainty, all eight evolutionary aims and the
large-data public release remain unfinished. Original full weighted timing
continues with sixteen workers and 1,131 checkpoints at 17:56 UTC; full weighted
biological fits remain zero and gated on its original completion/readback.
