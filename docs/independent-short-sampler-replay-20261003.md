# Full independent short-sampler output replay

The [observation-aware version-2 reader](independent-short-sampler-replay-v2-20261003.md)
has now passed full-grid software qualification and is separately queued for
all 1,620 original roles. It retains distinct unknown-residue draws and checks
both logs against concrete observations without asserting joint identity.
Version 1 remains frozen and its original jobs are unchanged; the historical
protocol below records its original strict cross-log contract. Full production
replay/readback/closure for both versions remain pending.

The complete 1,620-role independent output replay is queued behind the
original short sampler's full source/artifact/journal closure. It will decode
every successful role's three saved alignments and available site-property
records, preserve every failed or invalid original disposition, and verify
fresh residue-projection arrays by complete serialized readback. It launches
no inference. The full production replay has not started; ancestral posterior
adequacy and all eight biological aims remain incomplete.

## What the native logs contain

The generated model uses the same ancStates binding for ancestral alignment
output and category/state logging. Its category logger calls labeledNodeMap
on a tree whose internal labels have been dropped. Consequently the saved
FASTA contains ancestral residues, while the category JSON contains tips
only. The installed Bio/Alignment.hs defines that labeled-node restriction;
the generated programs and installed module hashes are bound by the
[immutable full plan](../metadata/independent_short_sampler_replay_plan_20261003.json).

Internal category labels are **unavailable** in this horizon. The replay does
not invent them from tip categories or claim a complete ancestral category
trajectory. A future full internal-category logger needs separate source and
native-software qualification before a new disjoint-seed inference horizon;
the running sampler and its frozen logger remain unchanged.

A pure native API probe verified the actual amino-acid state-index order
ARNDCQEGHILKMFPSTWYV. The projection alphabet is separately
ACDEFGHIKLMNPQRSTVWYX-, including unknown X and gaps. Both orderings remain
explicit; integer labels are not interchangeable or treated as scalar
evolutionary measurements.

## Full scope and admission

The original short sampler must first close every role, all 405 quartets,
135 effective inputs and 324 original aliases. Its completed archive, all
bound source/artifact hashes and complete disposition export are checked.
Failures do not prevent accounting for the full grid and never select a
smaller successful subset. The dependency wrapper verifies the original
sampler-closure handle and actual invocation-linked terminal journal before
starting this producer.

The [producer and reader](../scripts/run_independent_short_sampler_replay.py)
load the same immutable original job manifest, reconcile all chain IDs and
fresh seeds, and retain original configuration, command, process and receipt
bindings. The [separate replay engine](../scripts/independent_short_sampler_outputs.py)
uses the existing manual FASTA/Newick/projection implementations, importing
neither Biopython nor the original alignment block parser or anchor builder.

For every successful role it checks:

- The entire source/runtime rooted-clade correspondence and nonroot branch
  lengths at the unchanged strict tolerance, deriving all four candidate-node
  identities directly from descendant sets and the mapping table.
- Exactly iterations 0, 10 and 20, the full runtime node set, equal-width valid
  sequences, unchanged concrete observed residues, and the original input-X
  wildcard semantics. Every candidate length and mapping matches the original
  per-frame audit.
- All 21 scalar log rows and finite prior/likelihood/posterior arithmetic at
  the original tolerance. This is log-integrity arithmetic, not an independent
  likelihood computation or model-adequacy test.
- Every tip category/state pair: integer support, residue-order lengths and
  exact state-to-letter agreement with the corresponding ungapped native
  sequence. Missing tips, invented ancestors, duplicate JSON keys, extra
  fields and missing/duplicate/extra iterations fail.
- Every 4-by-20 rate-property cell: finite nonnegative values, category rows
  constant across amino-acid states, equal-weight mean rate one at 1e-10
  absolute/relative tolerance, and the original empty-condition schema.
  These checks do not independently reproduce gamma discretization or
  establish that the likelihood model is adequate.

Nonzero exits and original invalid-output dispositions retain their native
receipts and source hashes, with no successful arrays fabricated for them.
The replay does not silently reclassify or retry those attempts.

## Array schema, saved records and complete readback

Each successful role produces an uncompressed NPZ with exactly three arrays:

| Array | Type and axes |
| --- | --- |
| states | uint8; saved iteration × sorted source candidate node × concatenated extant-residue anchor |
| iterations | int32; exactly [0, 10, 20] |
| unanchored_residue_counts | int32; saved iteration × sorted source candidate node |

Anchor coordinates concatenate tips in sorted tip-label order, retaining
every nongap input residue in its within-tip order. Candidate nodes use
sorted source-node labels from the recorded candidate_mapping. The original
job manifest provides each input alignment and checksum. Homologous residues
anchored to several tips remain repeated coordinates, not independent
biological sites. Ancestral residues in columns with no extant residue are
counted separately; their identities are still preserved in the native FASTA.

The [versioned data dictionary](../metadata/independent_short_sampler_replay_data_dictionary_20261003.tsv)
defines role/frame fields, both state alphabets, array axes, coordinate order
and unavailable-category flags. It describes the schema, not completed
production outputs or biological acceptance.

Whole-role JSON checkpoints retain input/seed/alias identities, native receipt
links, all rooted-clade rows, candidate mappings, full per-frame alignment and
state digests, dimensions, lengths, unanchored counts, site-property digests
and explicit false scientific-eligibility/category-availability flags.
The full disposition export retains every requested role. All large native
and fresh-array artifacts remain outside Git, with paths/checksums in the
plans and eventual full archive.

The reader re-decodes all native samples and projections, compares every
serialized role/frame and every NPZ value/shape/dtype, and reconciles complete
scope totals and artifact hashes. Extra or missing role artifacts fail.
Completed producer/reader restarts are refused; incomplete arrays require
review. Whole-role checkpoints can be verified against fresh replay without
overwriting them. The reader shares the new decoder and is not a third
independent algorithm. Final closure requires all source/artifact hashes and
both actual original producer/reader completion journals.

## Software checks and preserved development failures

The [passed v3 software contracts](../metadata/independent_short_sampler_replay_software_validation_20261003_v3.json)
exercise the complete producer and reader with all 1,620 real role identities
and explicitly mocked source admission/native results. Two artificial failures
and all six intact roles in their unresolved quartets remain retained. Ten
restart/omission/identity/claim alterations and fifteen malformed native
property cases were rejected. These artificial outputs are software evidence,
not native production measurements.

Three existing completed native roles, one per prior, were separately decoded
read-only: nine actual saved alignments and 36 candidate frames. Every state
and unanchored count matched the original independent anchor-builder path.
A separate synthetic gap/X projection fixture exercises all 22 state letters.
One pure native alphabet API probe performed no alignment sampling or MCMC.
These checks precede the separate complete production replay; they do not
substitute for it or constitute a separate biological pilot.

Two unsuccessful software invocations remain preserved with their sources
and original journals. The first [supplied unequal-width synthetic sequences](../metadata/independent_short_sampler_replay_development_failure_20261003.json).
The second [confused original and independent statuses in its mocked loader](../metadata/independent_short_sampler_replay_development_failure_20261003_v2.json),
so the mock reported success without writing the expected array. Both
corrections concern software fixtures; native inference was unaffected.
The corrected v3 invocation has [actual terminal and source-hash evidence](../metadata/independent_short_sampler_replay_software_execution_20261003.json).

## Resources and execution

Each producer, reader and closer has two CPU equivalents, 16 GiB RAM, zero
swap and one BLAS thread. One role is decoded at a time. If every role passes,
the scope is 4,860 saved alignments, 19,440 candidate frames and 293,964,912
uint8 projected-state bytes. The largest projected array would contain
3,419,904 bytes. These are exact potential projection counts from the full
manifest; they do not bound JSON, FASTA or transient parser workspace.

Planning allows 8 GiB of output and 0.05–8 active hours per stage, excluding
the original sampler wait; neither is a calibrated guarantee or posterior
ETA. The free-disk guard requires 64 GiB and installed cgroups enforce RAM,
CPU and swap separately. At preparation, 306.87 GiB RAM and 10,375.01 GiB
disk were available. No GPU or paid resources are used.

At 07:48 UTC (03:48 EDT), [all six original replay/sampler handles and live caps were verified](../metadata/independent_short_sampler_replay_execution_checkpoint_20261003.json),
along with all 1,006 frozen pins. Replay producer/reader/closer controllers
2880769/2880773/2880777 were queued; zero replay roles had started. The
original sampler had 188 successful unclosed role checks. Its separate
resource observer was live with [fresh exact-handle evidence](../metadata/baliphy_sampler_resource_observation_execution_checkpoint_20261003_v4.json).
Historical allocation failures, full sampler/observer/replay closure,
longer adequate ensembles and biological root/model/predictor/calibration
qualification remain unresolved. All eight aims remain incomplete.

```bash
python scripts/record_independent_short_sampler_replay_checkpoint.py \
  --output NEW_SHORT_REPLAY_CHECKPOINT.json
```

Use new checkpoint names and inspect the original jobs; do not relaunch the
frozen producer/reader/closer. The [CPU environment](../environments/independent-short-sampler-replay-20261003.yml)
records dependencies. Run software checks only into a new audit directory and
receipt using scripts/check_independent_short_sampler_replay.py; this pure
native probe and read-only early QC remain distinct from production closure.

## October 3: ambiguity and joint logger correction

The [full input census](../metadata/independent_short_sampler_full_input_ambiguity_census_20261003.json)
found three `X` residues in one of 39 distinct alignment files, affecting 24
roles. A capped native software test shows separate FASTA/property loggers
can produce different valid conditional residue draws at an unknown observed
position. The frozen v1 decoder's exact cross-log tip equality can therefore
reject a distinct conditional draw. None of those roles had checkpointed at
the census time; this is not a demonstrated failure of their native jobs.
The full production replay still needs a separately qualified correction
that checks concrete observations, retains unknown-residue differences and
does not combine independent log representations. All roles remain required.

For future sources, the [joint logger correction](baliphy-joint-node-logger-20261003.md)
has passed all-prior native software checks with strict full-node state/letter
equality inside one record. All 405 future programs/1,620 disjoint seed roles
are prepared, but unlaunched. This cannot recover omitted internal category
labels or force joint identity between historical separate logs. The current
sampler, reader sources, limits and schedules remain unchanged.
