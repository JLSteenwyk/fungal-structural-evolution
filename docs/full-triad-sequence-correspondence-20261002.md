# Full sequence-derived duplicate/reference correspondence

The full sequence-only correspondence control is running for all 27,056
source-ready physical duplicate/reference triples. Each set is aligned using
both installed MAFFT and FAMSA under all six input permutations: **324,672
native alignment states**. This will support comparisons of measured structural
divergence when residue correspondence is determined from sequences, alongside
the completed structural-alignment correspondence alternatives. It does not
establish true homology, independent predictions or a duplication effect.

## Complete source scope

The independently closed sequence catalog retains all 31,235 original ordered
physical triples. All 27,056 source-ready triples are scheduled; the other
4,179 retain their original fields and explicit unscheduled disposition. The
catalog remains source-bound to all 283,409 original target contexts, 214,461
reference ties and 428,922 logical duplicate/reference sides. Missing models,
excluded parents, identical models and absent measurement-design edges remain
in the original source records. No favorable structural result selects inputs.
Source readiness is a work/availability condition, not biological duplication
or orthology qualification.

Every sequence is the full original protein string from the audited model/input
registry. Confidence filtering occurs after sequence alignment, preserving
original sequence offsets. Stable model indices are assigned by sorting the
three versioned model identities; every original A/B/reference role permutation
remains linked explicitly. Identical physical sets may be reused without losing
logical occurrences. The actual census found 27,056 distinct sets, a maximum
single-protein length of 2,271 and a maximum summed triple length of 3,819.
All actual source letters are among the 20 standard amino acids.

Full producer and separate original-source/role/sequence readback passed all
31,235 links and 27,056 sequence sets. Closure verified **321,761 source/artifact
hashes and both original producer/reader completion journals**.
[Completed catalog locator](../metadata/full_triad_sequence_catalog_completed_20261002.json).
The full source hash archive, sequences and original context records remain
outside Git. The software contract also checked physical reuse, role swapping
and unscheduled originals, rejecting ten rehashed altered exports.

## Native methods and exact output requirements

MAFFT 7.525 uses `--amino --anysymbol --thread 1 --auto`. The author's
[manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html) documents
automatic strategy selection; its
[symbol documentation](https://mafft.cbrc.jp/alignment/software/anysymbol.html)
explains unknown-character scoring under `--anysymbol`. FAMSA
2.5.2-2598410 uses `-t 1` with its default single-input MSA settings;
[author repository and usage](https://github.com/refresh-bio/FAMSA).
The actual executable, MAFFT helper files, versions/help and tested commands
are pinned in the [environment record](../metadata/full_triad_sequence_alignment_environment_20261002.json).
No package update or new structure prediction was made.

Both methods are run for all six permutations, without selecting an input order
or alignment based on a structural distance. Fully observed MSA columns map to
three original protein positions. Every ungapped output must exactly reproduce
its original full sequence, with all three unique identifiers, equal column
counts and no all-gap column. Numerical coordinate fitting is a separate stage.

The SQLite checkpoint preserves every job's exact input FASTA, raw output/error
bytes, checksums, intended native command, PID/create time, observed command when
available, elapsed time and return status. Fast native processes may exit before
their command can be sampled; that observation remains explicit. Timeouts,
nonzero exits and invalid alignments retain separate dispositions and unavailable
correspondence metrics. A 600-second native guard retains timeout records rather
than deleting those states.

The separate reader uses Bio.SeqIO and independent per-sequence position counters
to reconstruct every raw MSA and common triple. It verifies the complete
method/order/set grid, exact source sequences/roles, command/input/output hashes,
status and every exported position. Completion additionally requires full hashes
and both original producer/reader completion journals. Valid alignment output
does not verify optimality, alignment accuracy or correct homology.

Actual installed-tool checks passed all 48 software states: four indel, repeat,
identical and nonstandard-letter fixtures, both methods and six orders. All
15 changed checkpoint/grid exports were rejected. Exact replay of an interrupted
task's owned input also passed. These are software contracts within the full
workflow, not a biological pilot or proof that production is complete.
[Fixture evidence](../metadata/full_triad_sequence_alignment_fixture_validation_20261002.json).
A separate [all-gap-column contract](../metadata/full_triad_sequence_alignment_all_gap_contracts_20261002.json)
prepended one all-gap column to a valid native alignment while preserving all
ungapped sequences, identifiers and equal widths. Both independent parsers
rejected this isolated defect; original native outputs remain unchanged.

## Resources, progress and reproduction

Resource estimates preceded launch: two CPU workers, one native/BLAS thread
each, 32 GiB RAM, no swap, a 32 GiB output allowance and 100 GiB free-disk
reserve. A 64 KiB-per-state raw storage scenario is about 19.82 GiB before
compression and SQLite overhead. An illustrative 0.1–6 CPU seconds per state
gives a 4.5–270 hour active-time planning range with two workers; it is
uncalibrated, excludes source hashing/extreme cases and is not an ETA. No GPU,
paid resources or existing-job changes are used.
[Resource estimate](../metadata/full_triad_sequence_alignment_resources_20261002.json).

The [native runtime checkpoint](../metadata/full_triad_sequence_alignment_runtime_checkpoint_20261002.json)
verified the three original live handles and recorded 2,048/324,672 states,
all valid at that observation. Later progress stays in the run's atomic
`state.json`. The [project checkpoint](../metadata/project_runtime_checkpoint_20261002_v7.json)
rechecked 323,636 distinct closed bindings, 25 live pipeline handles, six
original scientific/retrieval jobs, 29 terminal successes and ten preserved
failures; background comparisons had reached 196,544/298,848 states.

```bash
python scripts/prepare_full_triad_sequence_catalog.py \
  --plan metadata/full_triad_sequence_catalog_plan_20261002.json
python scripts/readback_full_triad_sequence_catalog.py \
  --plan metadata/full_triad_sequence_catalog_plan_20261002.json \
  --output results/structural_comparisons/full-triad-sequence-catalog-20261002-v1/readback.json
python scripts/run_full_triad_sequence_alignments.py \
  --plan metadata/full_triad_sequence_alignment_plan_20261002.json
python scripts/readback_full_triad_sequence_alignments.py \
  --plan metadata/full_triad_sequence_alignment_plan_20261002.json \
  --output results/structural_comparisons/full-triad-sequence-alignments-20261002-v1/readback.json
python scripts/close_full_triad_sequence_stage.py \
  --plan metadata/full_triad_sequence_alignment_completion_plan_20261002.json
```

These commands identify the frozen completed/running/queued jobs. Do not start
duplicates. Reproduction uses new plans/output locations; restarting the native
run after its original process is authoritatively stopped reuses only checkpoints
whose configuration, raw hashes, source sequence and output fields match exactly.
The future native completion locator is
`metadata/full_triad_sequence_alignment_completed_20261002.json`.

## Required structural and evolutionary integration

The complete sequence-derived geometry workflow is now implemented and queued
behind the original native alignment/reader closure. It covers all **649,344
fit dispositions**: 27,056 ordered triples, two alignment methods, six input
orders and two confidence masks. Production preflight and geometric results
are still pending; the software checks below do not establish their completion.

## Full sequence-derived geometry queued October 2 UTC

A full producer and separate raw-MSA reader first count every native state,
masked core, unique eligible core, residue occurrence and required full-PDB
input/byte. Both original preflight journals and complete source hashes must
close before any coordinate fitting. This measures actual work across the
entire dataset, without a separate biological pilot.

Fits project sorted sequence indices back to the original A/B/reference roles.
Every pair uses the same original three-way residue matches. Confidence masks
use membership in the previously verified original-position lists; rounded PDB
confidence never determines retention. For example, a source confidence of
69.996 can print as 70.00 in a materialized PDB while remaining excluded from
the authoritative pLDDT70 mask. Outputs therefore distinguish
`joint_authoritative_plddt70_fraction` from `joint_plddt70_fraction` calculated
from rounded PDB values. Mean confidence values also reflect materialized PDB
precision. Filtering applies jointly across all three proteins after full
sequence alignment. Coverage uses each full original protein length; retained
mask length and its separate coverage are exported explicitly.

The producer computes proper SVD rotations for AB, AR and BR, signed AR−BR,
sequence identities, confidence summaries and rotation-uniqueness dispositions.
All six original rational coverage screens remain explicit, separately for
the new common core, inherited both-order three-edge eligibility, and their
intersection. The actual full 432,896-view structural mapping was checked:
each original triple/mask has constant both-order three-edge screening gates
across all eight structural-order views.
[Full gate-constancy evidence](../metadata/full_triad_sequence_fit_baseline_gates_verified_20261002.json).
This consistency does not promote excluded logical contexts or parents.
Native failures, short cores, source-rejected masks and degenerate geometry
remain rows with explicit status; uncomputed numerical metrics remain blank.

Transactional SQLite checkpoints commit complete triads. Restart reconstructs
every saved fit from its original immutable sources, checks compressed payload
hashes and exact calculated fields, then rebuilds the complete export without
duplicating rows. An exclusive run lock prevents concurrent writers. A completed
receipt refuses restart. The independent reader rebuilds original positions
from raw MSA using Bio.SeqIO, projects authoritative masks separately and uses
a quaternion eigenproblem to calculate every fit. It checks all numerical
fields at 1e−9 absolute/relative tolerance, RMSD metric bounds, every source/role/
mask/screen field, complete row ordering and the exact SQLite/export agreement.
Closure requires both original producer/reader completion/resource journals
and the complete source and artifact hash inventory.

Software contracts passed the full 72-state/144-fit grid across five physical
model sets and six ordered role triples, including physical reuse and a
reversed duplicate-role projection. Tests exercise real raw-MSA parsing,
SQLite and PDB reading, unequal lengths, nonconsecutive confidence masks,
the confidence-rounding boundary, native failure, a two-residue common core,
source rejection, collinear and reflected coordinates, inherited exclusions
and sign reversal. The independent maximum RMSD/contrast difference was
1.849e−14 Å. All 16 rehashed false geometry exports and three altered preflight
summaries were rejected. Recovery from an owned synthetic interruption rechecked
and reused 24 committed rows. Production source-closure I/O is stubbed only in
these software fixtures; production requires the real complete source closures.
[Software contract locator](../metadata/full_triad_sequence_geometry_fixture_validation_20261002.json).

Resource estimates preceded launch: two-CPU quota, 32 GiB RAM, no swap, one
BLAS thread and no GPU or additional charges. Geometry allows 32 GiB output
and requires 100 GiB free disk. Uncalibrated planning ranges are 0.25–8 hours
for each preflight/reader stage and 1–24 hours for each geometric fit/readback;
they exclude dependency waits and are not completion ETAs. Exact unique work
and input-byte demand must be verified by full preflight first.
[Resource estimate](../metadata/full_triad_sequence_geometry_resources_20261002.json).

The six original queued launch identities, creation times, exact commands and
actual cgroup limits are recorded in the
[launch inventory](../metadata/full_triad_sequence_geometry_followup_launches_20261002.json).
All waits use original invocation-linked journals, not collected unit defaults.
The [frozen workflow](../metadata/full_triad_sequence_geometry_followup_workflow_20261002.json)
and [launcher](../scripts/launch_full_triad_sequence_geometry_followup.py) describe
the queue; **do not launch duplicates**. The two source plans are
[preflight](../metadata/full_triad_sequence_fit_preflight_plan_20261002.json) and
[geometry](../metadata/full_triad_sequence_geometry_plan_20261002.json).
Future completion locators are
`metadata/full_triad_sequence_fit_preflight_completed_20261002.json` and
`metadata/full_triad_sequence_geometry_completed_20261002.json`.

## Remaining evolutionary integration

The complete comparison with structural correspondence alternatives is now
implemented and queued, as described below. Its full production results remain
pending. Every original logical context and parent exclusion must still be
integrated without promotion based on a favorable common-core result.

## Full correspondence comparison queued October 2 UTC

The new stage consumes only the complete, independently closed 649,344
sequence-derived and 865,792 structure-derived fits. Both original source
producer/reader completion journals and full source/artifact hashes are required.
Original sequence sets, A/B/reference roles, masks, structural mapping definitions,
native payload hashes and complete row grids are checked before comparison.
The same original ordered triple is compared throughout; no favorable alignment
order, reference, mask or core is selected.

| Output | Full planned row count | Scope |
| --- | ---: | --- |
| Sequence order robustness | 108,224 | 27,056 triples × two masks × two sequence methods; all six input orders retained |
| Correspondence comparison groups | 216,448 | Sequence groups × two structural core definitions |
| Paired comparison states | 10,389,504 | Every comparison group × six sequence orders × eight structural orders |

Each paired state retains both fit statuses, exact residue-triple hashes and
sizes, intersection and union counts, Jaccard overlap, correspondence equality,
numerical uniqueness and strict numerical contrast-sign agreement. Empty-union
Jaccard values remain unavailable. Every one of the 18 shared numerical fields
has its sequence-minus-structural difference exported where both values exist.
This includes residue count, original protein-length coverage, AB/AR/BR RMSDs,
AR−BR, sequence identities, rounded-PDB confidence summaries and rotation
curvatures. These differences describe sensitivity to correspondence and coverage;
the two cores can contain different positions and do not define independent
evolutionary effects. Uncomputed values remain blank. Six screens retain each
source's combined pass and their intersection for every order pair.

The separate sequence-order summary retains all six statuses, native statuses/
payload hashes, correspondence hashes, 19 numerical ranges (including the
authoritative source-mask fraction), all six core/inherited/combined screening
bitmaps and source exclusions. All-order qualification requires all six orders;
any-order qualification is reported separately. Joint groups retain the full
48-pair layout, numerical difference ranges, missing-value counts, correspondence-equality
counts, strict numerical sign agreement and all/any-pair screen bitmaps. A fully
passing group requires all 48 order pairs. Numerical signs use the exported
values and have no significance or evolutionary-polarity interpretation.

The producer streams complete triples, uses set intersections for correspondence
overlap and calculates full order-pair differences. A separate reader reconstructs
every sequence correspondence from raw MSA with Bio.SeqIO, calculates overlap
using a sorted-tuple merge and uses SQL for pair subtraction, ranges, means,
status counts and order/screen bitmaps. It checks every output field and complete
ordering, all missing values, every summary and the checkpoint/export agreement
at 1e−9 absolute numerical tolerance. No producer summary or comparison algorithm
is imported. Full source hashes and both original completion/resource journals
must close after the independent reader.

Deterministic compressed checkpoints retain each triple's complete sequence,
pair and joint-group exports. Interrupted restart regenerates and exactly checks
every existing checkpoint before rebuilding the full tables. An exclusive lock
is held through receipt creation. Complete receipts refuse restart; live jobs
must never be duplicated. The full software test used the existing actual raw-MSA/
PDB/SVD fixture sources and the entire six-by-eight grid across six ordered
triples, including role reversal, two masks, two methods, both structural cores,
incomplete/source-rejected/degenerate states and inherited exclusions. All 2,304
pair states/48 joint groups/24 sequence groups passed; 11 rehashed false pair
exports and six false group exports were rejected. Interrupted recovery checked
two committed triple blocks. The maximum independent numerical difference was
4.441e−16. Production closure I/O is stubbed in these software contracts;
the production workflow requires the real full closures.
[Software evidence](../metadata/full_triad_correspondence_comparison_fixture_validation_20261002.json).

Resources were estimated before launch: two-CPU quota, 32 GiB RAM, no swap,
one BLAS thread, 32 GiB output allowance and a 100 GiB free-disk reserve.
The uncalibrated planning range is 1–36 hours for each producer/reader, excluding
dependency waits; it is not a completion ETA. No GPU or new charges.
[Resource estimate](../metadata/full_triad_correspondence_comparison_resources_20261002.json).
The three original jobs wait behind the original sequence-geometry closure;
[launch identities](../metadata/full_triad_correspondence_comparison_launches_20261002.json),
[frozen plan](../metadata/full_triad_correspondence_comparison_plan_20261002.json)
and [workflow/launcher](../scripts/launch_full_triad_correspondence_comparison.py).
The future completion locator is
`metadata/full_triad_correspondence_comparison_completed_20261002.json`.
Do not launch duplicates. Independent reproduction requires new plans/output
locations or a checked restart after the original process is authoritatively stopped.

Full production native alignment, sequence geometry and correspondence comparison
are still pending. Original source/context/parent eligibility, joint-mask/method/
core robustness, fixed matched controls, domain/orientation/PAE and predictor
differences must still be integrated into qualified phylogenetic and calibrated
biological analyses. All 31,235 originals and 283,409 contexts/214,461 ties/428,922
sides remain linked through the closed source catalog, including unscheduled
and excluded originals.

The two aligners and six orders are dependent sensitivity alternatives, not
additional biological replicates. Sequence-derived correspondence reduces one
source of structural alignment dependence; predicted coordinates still derive
from sequences. Domain/orientation/PAE and prediction-source controls, biological
orthology/phylogenetic qualification, matched backgrounds and calibrated
uncertainty/multiple testing remain required. Signed A/B reference contrasts
do not establish ancestral polarity or which duplicate evolved faster.
All eight scientific aims remain incomplete; GPU prediction remains paused.
