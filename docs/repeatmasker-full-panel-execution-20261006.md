# Full-panel RepeatMasker execution

This is the required all-instance repeat annotation after complete de-novo
library discovery. It will submit all 526 taxa only after the independent
RepeatModeler audit verifies one classified consensus library per manifest
taxon. It does not use a discovery subset or select a smaller taxon panel.

For each taxon, a source-bound worker will verify the completed RepeatModeler
library and combine it with the fungi-and-ancestor curated Dfam consensus
extract. The latter contains nine Dfam 4.0 records in the current release; its
checksum and original FamDB audit are pinned in the configuration. Distinct
`RM2__` and `DFAM40__` prefixes preserve the source of every library record.

RepeatMasker 4.2.4 will use RMBlast with eight threads, sensitive mode,
lowercase masking, GFF, alignment, and no bacterial insertion-sequence check.
The installed FamDB version reads its data-directory setting from local
`famdb.conf`, not an environment variable. The versioned configuration helper
sets that path to the audited Dfam directory and writes a checksum receipt
([`famdb_local_configuration_20261006_v1.json`](../metadata/famdb_local_configuration_20261006_v1.json))
before RepeatMasker is eligible to launch.
Each output will retain the original source FASTA linkage, combined-library
hash, full command, `.out`, `.tbl`, `.gff`, and alignment outputs. Twelve
concurrent tasks reserve 96 CPU threads and 384 GB RAM on the 192-core CPU
capacity partition.

The annotations will quantify repeat classes and gene proximity only after
independent output readback and explicit assembly/annotation sensitivity
handling. Repeat annotation alone does not establish repeat-mediated genome
rearrangement or structural-protein evolution.
