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
