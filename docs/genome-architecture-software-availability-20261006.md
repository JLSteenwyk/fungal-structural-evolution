# Genome-architecture software availability

The original [software availability receipt](../metadata/genome_architecture_software_availability_20261006_v1.json)
records the pre-installation state. It is retained as historical evidence; the
isolated environments below now provide the required repeat and collinearity
software. No full-panel repeat or synteny analysis has started, so this is still
a launch-preparation record and not a biological result.

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

The repeat environment also solved and installed successfully in the project
cache on October 6. Its explicit, 206-package Linux lock is
[`genome-architecture-repeat-linux-64.lock.txt`](../environments/genome-architecture-repeat-linux-64.lock.txt).
RepeatMasker reports version 4.2.4, RepeatModeler reports version 2.0.9, and
RMBlast is available. The [Dfam source audit](../metadata/dfam_fungal_library_20261006_v1.json)
records Dfam 4.0 (FamDB format 3.0.0), including curated consensus and the
root/Fungi uncurated-HMM sensitivity partitions. Its source files remain outside
Git and are checksum-recorded in that receipt. The
[full-panel DNA input readback](../metadata/full_panel_repeat_input_readback_20261006_v1.json)
also independently re-read the source-bound compressed FASTA for all 526 taxa.

Before the full, non-pilot 526-taxon run, the remaining gate is a defensible
all-taxon resource/restart plan and an execution controller that preserves each
taxon's provenance and failure disposition. EDTA remains an optional third
repeat-method sensitivity rather than the default fungal repeat call, pending a
fungal-appropriate configuration and resource plan.
