# Full-panel RepeatModeler execution

This stage launches de-novo repeat-library discovery for every one of the 526
source-verified genomes. It is an all-taxon array, not a pilot: task indices
0–525 map exactly to the immutable input manifest whose fresh readback covers
26,645,610,508 bases and 1,070,884 records.

The Slurm array reserves eight CPU threads and 32 GB RAM per taxon, with at
most twelve concurrent taxa (96 threads and 384 GB reserved). It targets the
192-core `gpu` partition for CPU capacity; no GPU is requested by this stage.
This leaves capacity for the currently active structural-atlas readers and keeps
the run within locally available, no-cost resources. The seven-day task limit bounds
individual large assemblies. Failed tasks retain their work directory. The
controller uses RepeatModeler's recovery function after more than one completed
round; because the software itself rejects early-round recovery, earlier
attempts are retained under an `interrupted_` directory and restarted from the
same checksum-verified source rather than silently restarted from another
source.

Before any tool invocation, the controller performs a fresh SHA-256 check of
the compressed publisher-bound FASTA, decompresses it, and checks its record
and base census against the immutable manifest. It runs BuildDatabase with
RMBlast then RepeatModeler 2.0.9 with eight threads and a taxon-derived,
documented seed. A completed receipt requires one non-empty classified library
and its SHA-256. RepeatModeler discovery is not a genome-wide repeat annotation
and therefore cannot by itself support a claim about repeat content, synteny or
structural evolution.

The worker pins both layers of FamDB configuration: `FAMDB_DIR` identifies the
FamDB program installation for RepeatModeler, and its local `famdb.conf` points
to the audited Dfam data directory. If discovery has completed but embedded
classification did not create a classified library, the worker runs
RepeatClassifier directly on the completed consensus under that same
configuration. This preserves discovery output and prevents a configuration
failure from being interpreted as biological absence of repeats.

[`audit_repeatmodeler_full_panel_v1.py`](../scripts/audit_repeatmodeler_full_panel_v1.py)
independently audits all 526 expected receipt locations and verifies the source
and classified-library SHA-256 values of every completed task. It reports
completed, incomplete, absent, or invalid-receipt dispositions without changing
worker output; this provides the full-scope completion gate before RepeatMasker
instance annotation begins.

The optional LTR structural module is not enabled here because its complete
dependency chain was not qualified in the pinned environment. It remains a
separate full-panel sensitivity stage, not a reason to exclude any taxon or to
substitute a smaller analysis. The subsequent RepeatMasker stage will combine
the per-taxon de-novo library with the checksum-recorded Dfam source library
and will write all-instance annotations, class summaries and provenance
receipts.
