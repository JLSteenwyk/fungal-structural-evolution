# Full duplicate/reference three-way work design

This stage preserves all 283,409 original duplication contexts, both reference
designs, 214,461 tied-reference records and 428,922 logical duplicate/reference
sides. All 121,490 original availability-side links remain in the source
records. It designs three-way comparisons; it does not establish shared
residue correspondence, biological orthology or a duplication effect.

The full producer and independent SQL reader completed. The
[completion record](../metadata/full_reference_triad_work_design_completed_20260930.json)
binds both exact original process journals and the complete source/artifact
hashes. The
[unit-labeled count table](tables/full_reference_triad_work_counts_20260930.tsv)
is published only after that closure. All eight scientific aims remain
incomplete; GPU prediction remains paused.

## Source preservation and model identity

`context_triad_design.jsonl.gz` retains the entire previously closed
`context_measurement_design.jsonl` record as `source_design`. It adds the
primary duplicate AB work state and a reference-triad link for every original
tied reference under each design. Missing reference models get a null triad
identifier and an explicit exclusion. Empty reference sets remain empty; their
context denominators remain in the receipt and count table. No parent, native
assignment, lexical gene, missing structure or source field is changed.

Every modeled reference defines an ordered triple of versioned models with
roles A, B and reference. `ordered_model_triads.jsonl` deduplicates these
triples for physical work while retaining their logical occurrence counts.
The A/B orientation follows the original gene identities rather than the
queue's table position. Model versions are integers. Deduplication retains
A/B orientation because signed asymmetry depends on those roles.

Three distinct versioned models and three distinct model IDs are separately
recorded. Different versions of one model do not establish independent
proteins. Identity exclusions remain explicit for AB, A–reference and
B–reference. Three distinct IDs are a source identity screen, not proof of
biological copy independence or independent prediction errors. Identical
sequences and other source/annotation controls remain in the original record
and require downstream assessment.

## All three edges and source-only work readiness

The original primary pair catalog supplies AB. The complete full-reference
pair catalog supplies A–reference and B–reference, including its original
additional/existing-primary work labels. Each edge carries its pair hash,
actual current endpoint order and desired role direction: A to B, reference
to A, reference to B. The desired-direction order describes endpoint
orientation; it does not select a favorable native alignment or replace its
original source-order provenance. Identical and outside-design comparisons
have null endpoint/order fields; no alignment or zero distance is invented.

A reference link is source-ready for correspondence only if the original
parent is eligible, the primary duplicate comparison is queued, all three
edges are in their designated complete measurement catalogs, and the triple
passes both identity screens. Source-guide and both-guide native assignment
sensitivities are recorded separately. They do not change source readiness,
lexical choice or the original tied-reference set. A physical triple measured
for another context cannot promote an excluded parent.

The closed full design has **31,235 unique ordered model triples**: 28,686
have three versioned models/IDs, 2,547 have two, and two have one. There are
**27,056 triples with at least one source-ready logical link**. These define
**432,896 potential correspondence states** across both confidence masks and
all eight AB/AR/BR input-order combinations. This is a future full grid, not
completed native residue mapping or accepted statistical observations.

Parent-eligible lexical source-ready references assigned to both duplicates
under both native guides number:

| Source guide | Availability design | Sequence-first design |
|---|---:|---:|
| Profile | 26,283 | 19,318 |
| MAFFT | 26,265 | 19,300 |

These count original lexical reference records. Counts overlap across guides
and designs and must not be pooled. They precede measured pair coverage,
common-core coverage, prediction-error controls and phylogenetic inference.
All missing/excluded records remain in the full export.

## Independent validation and execution

The independent reader shares source/proof I/O and summary field names only;
it imports no producer projection functions. Keyed SQL primary/reference
catalogs reconstruct every edge, actual endpoint, pair key, direction,
identity exclusion and contextual flag. It compares every full source field,
checks context semantic identities and ordinal grids, exhausts every physical
triple through its logical links, and independently aggregates all counts
using SQL. Producer and reader completion journals are captured using exact
PID, creation time and command identity. Default systemd success fields
after transient-unit collection are not sufficient completion evidence.

Full software fixtures passed 16 contexts/30 tied references/60 logical
sides. Ten rehashed false exports were rejected, including parent/identity
promotion, missing-reference removal, lexical replacement, swapped directions,
model-version promotion, altered occurrence counts, removed empty contexts
and changed pair hashes. Fixture proof/journal stubs are synthetic software
inputs, not real-data qualifications or a sampling pilot. Actual source
preflight rechecked all 40 bindings of the prior closed full context design.

Plans, source preflight, captured commands and resource estimates have the
`metadata/full_reference_triad_work_design_*_20260930.json` prefix. Data remain
outside Git at `results/orthology/full-reference-triad-work-design-20260930-v1`.
Each serial stage used two CPUs, 16 GiB RAM and no swap; resources were planned
before launch. Four GiB output allowance and 100 GiB disk reserve were
specified. The 0.1–4 hour envelope was uncalibrated planning, not a native
alignment ETA. No GPU inference, repeated optimization or paid resource was
launched. Preserve immutable plans, scripts and completed outputs.

## Remaining three-way analysis

Require closed native measurement unions and both-order/original-length
coverage for AB, AR and BR. Reconstruct every native residue map using its
actual source direction, masked input sequence and original position map;
retain unavailable, failed, numerically excluded and degenerate states. Use
all eight order combinations under both masks without choosing a favorable
direction. Intersect reference-to-A and reference-to-B maps and separately
require cycle consistency with AB. Fit the three structures on exactly the
same triples of original residues, with full-protein coverage denominators
and proper-rotation geometry checks. Sequence-locked correspondences, domains,
PAE/orientation, prediction uncertainty, supported phylogeny, annotation and
inferential calibration remain required. A signed distance difference alone
does not establish directional evolution or selection.
