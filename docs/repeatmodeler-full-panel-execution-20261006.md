# Full-panel RepeatModeler execution

This stage launches de-novo repeat-library discovery for every one of the 526
source-verified genomes. It is an all-taxon array, not a pilot: task indices
0–525 map exactly to the immutable input manifest whose fresh readback covers
26,645,610,508 bases and 1,070,884 records.

The Slurm array reserves eight CPU threads and 32 GB RAM per taxon, with at
most twelve concurrent taxa (96 threads and 384 GB reserved). This leaves
capacity for the currently active structural-atlas readers and keeps the run
within locally available, no-cost resources. The seven-day task limit bounds
individual large assemblies; failed tasks retain their work directory and are
recovered by the same controller rather than silently restarted from another
source.

Before any tool invocation, the controller performs a fresh SHA-256 check of
the compressed publisher-bound FASTA, decompresses it, and checks its record
and base census against the immutable manifest. It runs BuildDatabase with
RMBlast then RepeatModeler 2.0.9 with eight threads and a taxon-derived,
documented seed. A completed receipt requires one non-empty classified library
and its SHA-256. RepeatModeler discovery is not a genome-wide repeat annotation
and therefore cannot by itself support a claim about repeat content, synteny or
structural evolution.

The optional LTR structural module is not enabled here because its complete
dependency chain was not qualified in the pinned environment. It remains a
separate full-panel sensitivity stage, not a reason to exclude any taxon or to
substitute a smaller analysis. The subsequent RepeatMasker stage will combine
the per-taxon de-novo library with the checksum-recorded Dfam source library
and will write all-instance annotations, class summaries and provenance
receipts.
