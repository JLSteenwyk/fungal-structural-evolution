# Whole-proteome structural domain interval registry

The registry links every one of the frozen catalog's 1,290,278 AlphaFold models
to the independently audited full Pfam candidate-architecture database by exact
sequence hash. It retains the 1,319,513 original species/protein/model links,
checking each against the annotation protein table as well as the model table.
No prediction, structural clustering or evolutionary inference is performed
by this stage.

`structure_domains.sqlite` contains four related tables:

- `models`: structure identifier/version, sequence hash, length, path, raw-hit
  count and policy-agreement flag, including models without hits.
- `segments`: every distinct retained annotation hit across policies, including
  Pfam accession/version, type, clan, alignment bounds, envelope bounds and HMM
  coverage. Repeated hits remain separate through their source hit identifiers.
- `policy_hits`: retained hit membership for all four competition policies,
  conservative-architecture flags and domain-interval candidate flags.
- `protein_links`: every original taxon/protein association to its model.

Coordinates are **one-based inclusive source-protein positions**, with alignment
and envelope boundaries retained separately. Each interval must satisfy
`1 <= envelope_start <= alignment_start <= alignment_end <= envelope_end <= protein_length`.
Inconsistent same-hit fields across policies, repeated hits within a policy,
unknown model sequences, mismatched protein sequences and invalid bounds stop
the build. Database keys preserve distinct repeated occurrences.

Conservative architecture requires policy agreement and absence of rank ties,
unresolved overlaps, candidate nesting, alignment overlaps and retained partial
HMM hits below 0.70 coverage. A domain-interval candidate additionally requires
Pfam type `Domain`, HMM coverage >=0.70 and alignment span >=30 residues. These
are screening choices, not validated structural boundaries. All other types,
intervals and policy alternatives remain available for sensitivity analyses.
No-hit models remain explicitly unannotated; they are not domain-loss calls.

This registry does not extract coordinates or apply residue pLDDT/PAE criteria.
Those checks, independent registry readback, domain-level structural comparisons,
clustering and phylogenetic interpretation remain required. Exact sequence
identity supplies the intended residue numbering; it does not independently
validate every CIF coordinate here.

## Execution and resources

Script: `scripts/build_structure_domain_registry.py`.
Plan: `metadata/whole_proteome_structure_domain_registry_plan.json`.
Output: `results/domains/whole-proteome-structure-domain-registry-20260923-v1`.
Unit: `fungal-structure-domain-registry-20260923.service`.

The plan pins source data and code. Its catalog and architecture checksums were
also matched to the previously completed catalog and family-domain bridge
receipts before launch. Limits are one CPU, 16 GiB memory and no swap. Planning
allows 10 GiB output, at least 100 GiB free disk and an uncalibrated 0.5–12 hours;
no GPU or paid resources are used. The process was confirmed live in
`metadata/whole_proteome_structure_domain_registry_launch.json`.

Fixtures checked union-hit deduplication with preservation of four policy
memberships, conservative exclusion on policy disagreement and rejection of an
out-of-bounds envelope. A completed producer receipt will still be labeled
pending independent readback. At this checkpoint the build is running and no
registry result is claimed.

## Independent full readback queued

`readback_structure_domain_registry.py` now waits on the exact live registry
producer and requires its completed receipt, plan binding and database hash.
It reconstructs every model row from the catalog, every interval and policy
membership from the source annotation JSON, and every species/protein link
from both the catalog and original annotation protein table. It also checks
full table counts, foreign keys and database integrity. Source and output hashes
are checked before and after the readback.

The auditor does not import the producer's interval logic. In particular,
alignment overlap is recomputed by interval intersection, and partial-hit
exclusion is recomputed from each annotation's HMM coverage, instead of trusting
the producer's architecture summary counts. Fixtures show that stale summary
fields do not hide actual overlaps or partial HMM matches. This remains a check
against the same source annotations, not an independent Pfam search or a
validation of biological domain boundaries.

Plan: `metadata/whole_proteome_structure_domain_registry_readback_plan.json`.
Output: `results/domains/whole-proteome-structure-domain-registry-readback-20260923-v1`.
Unit: `fungal-structure-domain-readback-20260923.service`.
The controller is live; its launch record is archived in metadata. Limits are
one CPU, 16 GiB RAM and no swap; output planning is 0.01 GiB and uncalibrated
runtime planning is 0.5–12 hours after the producer. No GPU or paid service is
used. Independent validation is queued, not complete.

## Completed registry and independent readback

Both stages finished successfully. The full registry contains **1,290,278
models**, **1,319,276 distinct model/hit intervals**, **5,251,715 policy/hit
memberships** and all **1,319,513 protein links**. Of the models, 456,081 have
no qualifying Pfam hit; their records remain present.

The independent reconstruction matched every model field, interval, policy
membership, candidate flag and protein link, including recalculated interval
overlap and partial-HMM exclusions. Candidate Domain intervals number 594,425
under alignment/E-value ranking, 594,127 under alignment/bit-score ranking,
592,574 under envelope/E-value ranking and 592,279 under envelope/bit-score
ranking. These overlapping policy counts must not be added as unique domains.

The producer completed in about 264 seconds and the independent readback in
about 271 seconds, excluding their initial checks/waits. The SQLite SHA256 is
`90a2c795d2a0286fad037dad79a881757f365baaad50815ca6fb9fa4a4d0373e`.
Bound receipts are archived as
`metadata/whole_proteome_structure_domain_registry_completed_receipt.json` and
`metadata/whole_proteome_structure_domain_registry_completed_readback.json`.
The earlier queued/running text records launch history. Coordinate extraction,
residue confidence and PAE assessment, domain-level comparisons and evolutionary
interpretation remain unfinished.

## Complete domain extraction manifest

The union of candidate Domain hits across all four policies contains **594,797
model/hit pairs across 427,255 source models**. Preparing both alignment and
envelope spans yields **1,189,594 boundary associations**. Deduplication by the
full protein sequence hash and inclusive start/end coordinates retains
**1,078,592 unique intervals**, totaling **168,431,396 interval residues**; the
longest interval is 1,103 residues. The source models contain 198,306,563
full-protein residues.

The manifest retains every model/hit/boundary association in a deterministic
compressed table, so identical coordinates are stored once without losing
annotation provenance. It includes the full union, not a subsample. Interval
IDs hash the full source-sequence identity and bounds; they do not claim a
new evolutionary family or unique domain sequence.

Reproduce with `scripts/prepare_domain_extraction_manifest.py`, using a fresh
output directory. Outputs are under
`results/domains/domain-extraction-manifest-20260923-v1`; the archived receipt is
`metadata/domain_extraction_manifest_completed_receipt.json`. The producer
checked every boundary association against its original hit. The separate
`readback_domain_extraction_manifest.py` reconstructs the complete distinct
interval set with SQL UNION, then checks every exported model, source sequence,
path, bound, length and interval hash.

Coordinate extraction has not launched. Its dimension-based resource estimate
is recorded in `metadata/domain_coordinate_extraction_resource_estimate.json`:
four CPU workers, 32 GiB memory, no swap or GPU, 500 GiB output allowance and an
uncalibrated 6–96 hour planning range. The allowance covers potentially roughly
80–140 GiB of uncompressed all-atom PDB text, archives and validation outputs;
actual usage is not yet measured. Source coordinate hashes, full polymer
sequences, residue identity, missing atoms and confidence must be checked during
extraction. Domain boundaries, PAE qualification and downstream comparisons
remain separate checks.

The independent full manifest readback **passed** for all 1,078,592 intervals;
its receipt is `metadata/domain_extraction_manifest_readback.json`.

## Full coordinate extraction launched

The full 1,078,592-interval extraction is now running under
`fungal-domain-coordinate-extraction-20260923.service`; this supersedes the
not-yet-launched checkpoint above. The plan and verified live process record
are `metadata/domain_coordinate_extraction_plan.json` and
`metadata/domain_coordinate_extraction_launch.json`.

`extract_domain_coordinates.py` partitions all 427,255 source models into
1,000-model jobs and runs four CPU workers. Each worker reads each CIF once
into memory, checks that exact byte buffer against the catalog checksum,
verifies full polymer identity and atom-to-residue identities, and checks
finite coordinates, occupancy and confidence. Unsupported residues, alternate
locations, multiple models/chains, duplicate atoms and incomplete Cα coverage
receive explicit rejection dispositions. Missing other backbone atoms are
recorded on exported intervals and require downstream exclusion/review.

Exports retain all source atoms within the interval, with one-based local
residue numbering and the original start/end mapping in the manifest. PDB
coordinate/confidence rounding, serialized atom identity and fragment sequence
are checked before each export. Both boundary alternatives are retained through
the earlier reversible manifest. Each job produces an uncompressed tar archive,
a per-interval JSONL disposition/quality table and checksummed receipt. The
manifest includes source hashes, fragment sequence hashes, PDB hashes, mean
Cα pLDDT and fraction of Cα residues at pLDDT >=70. Those summaries do not make
an interval confidence- or PAE-qualified.

Independent Bio.PDB parsing checked serialized residue numbering, coordinates
and confidence on a real source fragment. Fixtures also rejected altered source
checksums and verified tar membership and completed-shard checksum checks. A
full independent archive readback is still required; these fixtures are not a
replacement for it. Completed shards carry receipts; unreceipted partial shards
require explicit recovery review rather than silent reuse. The top-level run
requires a fresh directory.

The run is capped at four CPUs, 32 GiB RAM and no swap. It requires 2 TiB free
disk before starting and checks a 1 TiB reserve before each source model. The
500 GiB storage allowance and 6–96 hour runtime range remain planning estimates.
On a worker failure, queued work is cancelled and already running bounded jobs
finish or reach their own checks. No GPU prediction, paid infrastructure,
biological domain-boundary validation or evolutionary inference is performed.

## Full domain archive readback queued

`readback_domain_coordinate_archives.py` now waits for the exact extraction
process and requires its completed receipt and plan binding. Before checking
archives it matches every job model to the original catalog and every job
interval to the full verified manifest, rejecting duplicate or missing scope.
It then uses four CPU workers to read every completed archive and disposition
record, checking hashes, exact tar membership and each exported atom against
original CIF arrays. Checks cover source/local residue numbering, atom and
residue names, elements, coordinates, occupancy, confidence, fragment sequence,
missing-backbone positions and confidence summaries. Standard PDB rounding
is allowed explicitly; altered coordinates are rejected by a fixture.

The auditor reconstructs expected atom sets and decodes PDB fixed-width fields
without importing the producer's extraction/serialization functions. The CIF
lexical parser is shared, so this is not a wholly independent file parser.
Rejected intervals must retain known identities and reasons; their rejection
causes are not independently adjudicated by this pass. The resulting receipt
will distinguish checked exports from excluded intervals. No PAE qualification
or biological boundary inference is implied.

Plan: `metadata/domain_coordinate_archive_readback_plan.json`.
Output: `results/domains/domain-coordinate-readback-20260923-v1`.
Unit: `fungal-domain-coordinate-readback-20260923.service`.
The live controller is recorded in `metadata/domain_coordinate_archive_readback_launch.json`.
Fixtures passed a complete real-source archive and rejected a changed coordinate.
Limits are four CPUs, 32 GiB memory, no swap and no GPU; output planning is
0.01 GiB and runtime planning is an uncalibrated 6–120 hours after extraction.
The readback is queued and not yet a completed validation result.
