# Joint ancestral sequence and category logging

The corrected logger has passed capped native software checks across all
three priors. The complete future source grid is prepared: **405 model
programs and 1,620 fresh seed roles**, covering all 135 effective inputs and
324 original configuration aliases. These sources have not been launched.
No posterior or biological aim is qualified by this software result.

## Why separate logs cannot be paired by iteration alone

The first all-node logger test failed strict sequence/state equality in five
records. Four involved ancestors and one involved an ambiguous observed
residue. Both native attempts exited zero; paired scalar, runtime-tree and
legacy FASTA files were identical. The original sources, outputs and failure
journal are retained in the [error recheck](../metadata/baliphy_full_node_logger_error_recheck_20261003.json).

The inspected primary source defines conditional ancestral sampling and
evaluates each registered logger through a separate action. The native tests
show that separately logged representations can contain different conditional
residue draws at the same iteration. This explanation is supported by the
co-evaluation comparison; it is not a trace of every runtime cache or random
number operation. [Pinned primary source evidence](../metadata/baliphy_joint_node_logger_primary_source_evidence_20261003.json).
The relevant upstream code is
[context logger execution](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/src/computation/context.cc)
and [ancestral sampling](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/haskell/SModel/Likelihood/VariableA.hs).

The future logger encodes the ancestral alignment in the **same record and
action** as the full-node category/state arrays. All nongap letters must match
their recorded native state indices exactly. The old separate alignment file
remains a preserved draw and must not be fused with the joint category record.

## Versioned changes and checks

Three exact, reversible edits apply to every existing model program:

1. Give the category logger the same internal-node label view used by FASTA.
2. Add `alignmentLines` to that category record, computed from the shared
   ancestral alignment expression.
3. Import `Data.Text` for the pure conversion into newline-free FASTA lines.

The inverse transformation recovers each original program byte-for-byte.
Observation distributions, likelihoods, indel/substitution models, priors,
reference initialization and transition definitions are unchanged. All
1,025 distinct full-grid input labels were checked for characters unsafe for
the installed JSON string encoder. This is a checked dataset restriction,
not a general JSON escaping fix.

Six sequential five-tip, 20-iteration native software runs compare original
and corrected sources at identical seeds under broad, centered and package
priors. These are computational checks, not a fungal pilot. The
[successful qualification](../metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json)
proves, for these runs:

- Exact paired bytes for all scalar TSV/JSON records, column maps,
  runtime trees and legacy FASTA files; unchanged rate/condition records.
- Exact joint sequence/state agreement for all nine saved frames, including
  36 ancestral-node records and 233 ancestral residue/category pairs.
- Matching gap patterns between legacy and joint alignments, with 18
  sequence differences explicitly retained across the two representations.
- Both original and joint tip draws satisfy every concrete observed residue.
  An input `X` leaves the residue uncertain; it does not permit a mismatch
  between states and sequences within a joint record.
- Four changed tip category records are retained. Exact equality of these
  conditional allocations across different logger evaluation scopes is not
  assumed. This does not establish their posterior distribution or mixing.
- Rejection of 33 altered records, including missing/invented ancestors,
  invalid categories, changed residue letters, lengths, iterations, missing
  alignments, unsafe lines and alignment/state disagreements.

The original software invocation completed with 57.582527 CPU seconds,
851,341,312 bytes peak memory and zero swap. Its actual invocation-linked
completion journal, source hashes and launch limits were verified;
collected systemd default status was not used as completion evidence.
[Execution proof](../metadata/baliphy_joint_node_logger_v4_software_execution_20261003.json).
The group had one CPU, 16 GiB RAM and zero swap. Each native run had 12 GiB
address space, 300-second CPU/wall bounds and a 64 MiB per-file cap. Six
sequential runs allowed at most 30 native wall minutes; planning allowed
1 GiB output and required 64 GiB free disk. No GPU or charges were incurred.

Two intermediate failures are preserved. Embedding multiline text directly
produced invalid native JSON because its newlines were not escaped
([v2 failure](../metadata/baliphy_joint_node_logger_v2_development_failure_20261003.json)).
The line-list version produced valid, internally consistent records, but its
paired test incorrectly required unchanged conditional category allocations
([v3 failure](../metadata/baliphy_joint_node_logger_v3_development_failure_20261003.json)).
The v4 test retains those differences and keeps strict same-record
sequence/state equality, observed residue constraints and complete paired
scalar/legacy alignment comparisons. Failed sources and native attempts were
not edited or restarted.

## Complete future grid and outstanding work

The [future source/seed inventory](../metadata/baliphy_joint_node_logger_future_models_20261003_v4.json)
contains all 405 programs, four roles per model, 540 roles per prior, and
1,620 distinct seeds disjoint from 3,243 earlier real/software seeds. It
preserves every original alias and source identity. The generated programs
total 3,568,815 bytes. It supplies neither execution commands nor a production
horizon. No current job source, native attempt, limit or schedule changed.

At 08:30 UTC the original short sampler had 747 successful unclosed role
checks; its queued independent replay had started zero roles.
[Verified original handles](../metadata/independent_short_sampler_replay_execution_checkpoint_20261003_v5.json).
The original resource observer remains live with
[fresh handle/limit evidence](../metadata/baliphy_sampler_resource_observation_execution_checkpoint_20261003_v5.json).
The broader covariance/timing queue was also
[reverified](../metadata/full_shared_entity_timing_execution_checkpoint_20261003_0831.json);
production fitting has not launched.

A full input census found three ambiguous residues in one of 39 distinct
alignment files, affecting 24 of the 1,620 original roles. None of those roles
had a checkpoint at the census time. The frozen queued replay v1 requires
exact equality of separately logged tip states and FASTA even at those
positions; it may reject distinct valid conditional draws. This is an open
reader issue, not an observed failure of those 24 native jobs. A separate
version must check concrete observed constraints, preserve disagreements at
unknown residues and avoid asserting joint identity between separate logs.
No affected role may be dropped or silently repaired.
[Complete ambiguity census](../metadata/independent_short_sampler_full_input_ambiguity_census_20261003.json).

Future source startup, full native-grid logger qualification, sampler/resource
closure, historical allocation failures, longer-chain memory/mixing,
independent likelihood checks and biological root/model/predictor uncertainty
remain open. Accepted phylogeny, reconciliation, dating, calibration and
**all eight biological aims remain incomplete**. GPU prediction remains paused.

## Reproduction

Use fresh output directories and receipt names. Preserve existing attempts:

```bash
python scripts/check_baliphy_joint_node_logger_v4.py \
  --output NEW_SOFTWARE_AUDIT --receipt NEW_SOFTWARE_RECEIPT.json
```

Run with the documented one-CPU/16-GiB/no-swap group and one-thread BLAS
settings. Source preparation also requires the verified original software
completion record; inspect
[the preparation script](../scripts/prepare_baliphy_joint_node_logger_v4_models.py)
before adapting destinations. Neither command authorizes posterior inference.
The [data dictionary](../metadata/baliphy_joint_node_logger_v4_data_dictionary_20261003.tsv)
specifies native alphabet order, coordinate interpretation and legacy/joint
draw distinctions. The Python environment remains
[the existing replay environment](../environments/independent-short-sampler-replay-20261003.yml).
