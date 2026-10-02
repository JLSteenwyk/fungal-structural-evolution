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

After native alignment/readback/journal closure, the next stage must fit all
three pair distances on the same sequence-derived positions under full and
pLDDT70 masks: 649,344 intended fit dispositions before within-triple coordinate
reuse. Confidence exclusions must apply jointly across all three proteins.
Proper-rigid geometry, original-protein coverage, all six existing screens,
source-parent exclusions, native failures and every original context must remain
explicit. These new fits must be compared with the existing structural
correspondence alternatives, including incomplete and conflicting cases.

The two aligners and six orders are dependent sensitivity alternatives, not
additional biological replicates. Sequence-derived correspondence reduces one
source of structural alignment dependence; predicted coordinates still derive
from sequences. Domain/orientation/PAE and prediction-source controls, biological
orthology/phylogenetic qualification, matched backgrounds and calibrated
uncertainty/multiple testing remain required. Signed A/B reference contrasts
do not establish ancestral polarity or which duplicate evolved faster.
All eight scientific aims remain incomplete; GPU prediction remains paused.
