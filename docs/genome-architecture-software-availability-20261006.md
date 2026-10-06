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

## Isolated environment specifications prepared October 6

The local Bioconda channel metadata available on October 6 includes the packages
needed for a two-track installation: RepeatMasker 4.2.4 and RepeatModeler 2.0.9
for de-novo/library-based repeat annotation, plus MCScanX 1.0.0 and JCVI 1.6.7
for independently implemented collinearity analyses. The corresponding pinned,
isolated environment specifications are
[`genome-architecture-repeat.yml`](../environments/genome-architecture-repeat.yml)
and
[`genome-architecture-synteny.yml`](../environments/genome-architecture-synteny.yml).

They deliberately do not combine all tools: repeat-library/database handling and
collinearity have different external inputs and resource profiles. The synteny
environment solved and installed successfully in the project cache on October 6;
its explicit, 351-package Linux lock is
[`genome-architecture-synteny-linux-64.lock.txt`](../environments/genome-architecture-synteny-linux-64.lock.txt).
The installed MCScanX binary reports its normal usage and JCVI imports as 1.6.7.
No genome has been analyzed by either tool.

The repeat environment has not been installed. Before the full, non-pilot
526-taxon run, the remaining gate is repeat-environment installation together
with repeat-library provenance, all-taxon resource bounds, restart behavior, and
an independent input readback. EDTA remains an optional third repeat-method
sensitivity rather than the default fungal repeat call, pending a
fungal-appropriate configuration and resource plan.
