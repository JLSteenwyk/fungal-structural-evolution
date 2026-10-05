# Full selected-assembly DNA acquisition — October 5, 2026

Assembly DNA is needed to check CDS/annotation correspondence, genomic flanks
and apparent duplicate genes, rather than relying on protein annotations alone.
It also supports contiguity/ambiguity and contamination-sensitive reviews. This
stage retains the full fixed **501 fungal entries plus 25 outgroups**; it does
not change sampling, resolve species identity or accept any duplication event.

## All selected sources and prelaunch estimates

The complete 526-row analysis manifest is joined to existing assembly/publisher
provenance. **519 NCBI taxa** have exact recorded genomic FASTA URLs, matching
their versioned protein/assembly directory: 499 fungi and 20 outgroups. All
519 deposited assembly-statistic reports are verified and pinned. Their
whole-assembly totals sum to **26,331,168,114 bases**; the largest is
1,488,883,982 bases. Primary-assembly totals are not substituted.

All **seven external genomic FASTAs** are already cached, including both fungi
and five outgroups. Each genome URL matches one unique original publisher
receipt/file. Figshare file names, publisher MD5s and exact article versions
are matched to the archived source metadata. These sources are reverified and
read, not downloaded again. Existing external measurements already cover all
seven; this new stage retains their provenance and per-record sequences in a
common whole-dataset inventory.

The complete-source plan binds 1,060 original sources. It estimates 5–25 GiB
of compressed NCBI genomes from 26.33 GB of uncompressed bases; compression
and network throughput are uncalibrated. The planning range is 2–72 hours,
with a seven-day hard CPU/wall cap, not a finish ETA. Acquisition uses two
CPUs/16 GiB/no swap, two HTTP workers, a 128 GiB cumulative payload allowance,
8 GiB per-genome file limit, 160 GiB output allowance, 2 TiB initial free-space
requirement and 1 TiB emergency reserve. Existing local resources are used;
there is no paid provisioning, new prediction or GPU use.

- [Complete immutable source/resource plan](../metadata/full_assembly_dna_plan_20261005_v1.json).
- [Observed headroom and execution resources](../metadata/full_assembly_dna_resources_20261005_v1.json).
- [Original native/API launch](../metadata/full_assembly_dna_launch_20261005_v1.json).

## Exact genome acquisition and original-record inventory

For every NCBI entry, retrieve the checksum listing and genomic FASTA from
the exact original assembly-version directory. Require one unambiguous
publisher filename/checksum, verify MD5, preserve SHA-256 and response length,
and retain checksum-listing bytes. At most three bounded HTTP attempts use
separate retained files; there is no automatic native-stage retry. Two
workers share a synchronized payload-byte budget and emergency disk checks.
Resource/budget failures stop the stage; source/format errors remain explicit
per-taxon dispositions. Existing original files are never overwritten.

Stream every deposited FASTA record through EOF, including gzip trailer/CRC,
requiring unique nonempty IDs and valid IUPAC DNA. Retain original sequence
case and ambiguity, complete descriptions, sequence and uppercase-sequence
hashes, and per-record lengths. No gaps are split and no organelle or other
record is removed. Report whole-FASTA base counts, lowercase masking,
N/other-ambiguity content and record N50/L50. These are **FASTA-record** metrics,
not a substitute for NCBI contig N50 or evidence of biological ploidy.

Compare recomputed total length with deposited `all` scope, keeping every
difference as an assembly-scope review disposition. Agreement establishes
length consistency, not correct genes, independent species identity or absence
of contamination/haplotigs. Publisher-version matching binds the assembly file;
existing taxonomic exceptions and annotation uncertainties remain unresolved.

Twelve offline literal controls pass exact DNA/case/hash/ambiguity/N50, malformed
IDs/empty/non-IUPAC/header ordering, truncated gzip, missing/duplicate checksum,
mocked full two-request receipt writing and resource-budget rejection. There
are zero live requests in those controls and no corpus pilot.

- [Acquisition implementation](../scripts/retrieve_full_assembly_dna_v1.py).
- [Literal contract controls](../metadata/full_assembly_dna_fixture_20261005_v1.json)
  and [actual original execution closure](../metadata/full_assembly_dna_fixture_transport_20261005_v1.json).

## Current evidence and required independent check

The original acquisition is complete across **526/526 entries**, with zero
source/format errors: 519 downloaded NCBI genomes and seven publisher genomes
checked in place. The retained FASTAs contain **26,645,610,508 bases**;
**8,101,966,631 bytes** were downloaded. All 519 recomputed NCBI lengths match
the original deposited whole-assembly reports. Original API session 63351
returned zero; the exact native wrapper, original invocation journal and
all declared source/output bindings close.
[Full acquisition receipt](../metadata/full_assembly_dna_20261005_v1.json)
and [original execution closure](../metadata/full_assembly_dna_transport_20261005_v1.json).

The complete independent reader launched only after that closure under original
API session **98777**, with four CPUs, 32 GiB RAM and no swap. At **22:26 UTC**, it
had independently reconstructed **423/526** genomes and remained verified live.
Full completion and original API/native/journal closure are still required.
[Original reader launch](../metadata/full_assembly_dna_readback_launch_20261005_v1.json)
and [exact runtime checkpoint](../metadata/full_assembly_dna_readback_checkpoint_20261005_goal_2227.json).

It uses a separate Biopython sequence parser and raw header/whitespace passes,
checks every contig digest, reconstructs record/base/case/ambiguity/N50 metrics
with integer/rational arithmetic, and rechecks every publisher listing/version/
receipt. Fractional metrics use the unchanged explicit two-ULP roundoff bound;
integer/count/hash checks are exact. Seven literal controls reject altered
counts/metrics/contig hashes. The planning range was 2–72 hours, uncalibrated,
with seven-day safety caps. No corpus pilot or source normalization occurs.

- [Independent full reader](../scripts/readback_full_assembly_dna_v1.py)
  and [literal controls](../metadata/full_assembly_dna_readback_fixture_20261005_v1.json).
- [Full reader plan and completion gate](../metadata/full_assembly_dna_readback_plan_20261005_v1.json)
  and [resources](../metadata/full_assembly_dna_readback_resources_20261005_v1.json).

Plans and execution configuration retain exact commands and SHA-256 source
bindings. Large genomes and per-contig/disposition files remain outside Git
under `results/assembly-dna-20261005-v1`; external files remain at their recorded
original paths. Existing output roots reject re-execution; new runs require
separate versioned plans and roots. Cached source files are not normalized
or replaced to make checks pass.

After full acquisition/independent closure, join genome IDs and coordinates to
the exact GFF/CDS/gene/protein mappings, check coding-sequence/translation
correspondence, then evaluate genomic-flank and sequence evidence for suspected
annotation/haplotig/duplicate artifacts. These downstream checks and lineage
completeness/contamination/taxon review are still required. Genome acquisition
does not infer duplication, selection, ecology or structural acceleration.
All eight evolutionary aims and final project deliverables remain unfinished.
