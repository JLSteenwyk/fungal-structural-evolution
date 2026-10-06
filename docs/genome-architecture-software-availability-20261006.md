# Genome-architecture software availability

The full-panel repeat and synteny stage is not launched because the local
[software availability receipt](../metadata/genome_architecture_software_availability_20261006_v1.json)
reports no installed RepeatMasker, RepeatModeler, EDTA, MCScanX, i-ADHoRe or
JCVI executable. This is a concrete launch dependency, not a biological result.

`minimap2`, BLAST+, DIAMOND and MMseqs2 are installed and version-pinned in the
receipt, but general sequence alignment/search does not replace repeat annotation
or orthology-anchored collinearity inference. The workflow will select a
versioned repeat/collinearity implementation, record its database provenance and
resource/restart plan, and confirm no paid infrastructure before a full 526-taxon
launch. No substitute or reduced-scope analysis has been started.
