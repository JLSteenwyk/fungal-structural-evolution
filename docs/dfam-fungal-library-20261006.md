# Dfam fungal repeat-library provenance

The repeat workflow uses Dfam 4.0 in FamDB format 3.0.0 as a consistent
library source. The source audit is
[`metadata/dfam_fungal_library_20261006_v1.json`](../metadata/dfam_fungal_library_20261006_v1.json).
It records the Dfam release date (2026-05-22), source URL, exact SHA-256 and
byte count for each downloaded HDF5 partition, and output from the official
FamDB fungal-partition check.

The downloaded material is intentionally limited to the Dfam root database,
curated consensus sequences, and root plus Fungi (taxon 4751) uncurated-HMM
partitions. Curated consensus supports the primary library-based annotation;
the HMM partitions are retained for a separate fungal sensitivity workflow.
The library files are stored outside Git under
`data/references/dfam-4.0-curated-consensus-20261006/`.

The official downloader validated each archive MD5 before decompression. The
project audit independently hashes the decompressed HDF5 files and verifies
that the Fungi query finds the expected root and fungal uncurated-HMM
partitions. This establishes retrieval and provenance only. It is not a repeat
annotation, and it does not establish repeat identity, repeat-mediated genome
change, synteny, or structural-evolution associations.
