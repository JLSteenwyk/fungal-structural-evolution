# Full candidate domain-architecture representation

**September 22 update:** construction is complete for all 5,713,603 queries and
5,815,847 protein links across 526 taxa. Independent full readback is running.
The representation contains 5,767,208 per-query retained-hit alternatives;
these are not globally distinct biological architectures. The producer receipt
is archived at `metadata/full_candidate_domain_architecture_completed_receipt.json`.
The launch descriptions below preserve the initial state and planning context.

The complete clan-competition output passed independent readback over all
5,713,603 queries and 8,103,610 annotation hits. There are 53,248 queries whose
retained hit sets differ across the four policies, 1,534,637 raw hits with
HMM coverage below 0.70, and 1,991,329 queries with no gathering-threshold hit.
These last queries are not confirmed biological domain absences. Production
and independent audit receipts are archived under
`metadata/full_domain_competition_completed_{receipt,readback}.json`.

The next full-scale stage began September 22. It joins every query's audited
competition outcome to raw annotation fields and creates a candidate
architecture database, with the original 5,815,847 protein links across all
526 taxa. Four marker-only queries also remain in the query universe. A
sequence shared by multiple proteins keeps every taxon-specific protein link.

The representation retains all Pfam annotation types, including Family and
Repeat; they are not silently relabeled as structural domains. Repeated model
occurrences retain multiplicity and coordinates. For each distinct retained
hit set, annotations are sorted by alignment start/end and deterministic hit
ID. The four policy records point to these alternatives and retain counts of
primary score ties, unresolved overlaps and candidate nesting relationships.
The complete per-hit suppression/blocker history remains in the audited source
competition table, identified by the source receipt.

Every alternative also records all alignment-overlap pairs, partial-HMM hit
counts, ordered model/type tokens and their signature. Coordinate sorting
does not resolve nesting or establish an unambiguous biological domain order.
Token signatures are descriptive annotation identifiers, not homology or
functional classifications. Empty annotations have an explicit no-hit state.
This layer is intended to support later family-aware domain comparisons, not
to supply domain gains/losses without further validation and phylogenetic work.

Implementation: `scripts/candidate_domain_architecture.py` and
`scripts/build_candidate_domain_architectures.py`. The deterministic fixtures
in `scripts/check_candidate_domain_architecture.py` verify repeated-model and
annotation-type preservation, partial and overlapping matches, alternative
policies, explicit missing hits and invalid-ID rejection.

The plan is `metadata/full_candidate_domain_architecture_plan.json`; the unit
is `fungal-candidate-architectures-20260922.service`. It uses one CPU equivalent,
8 GiB memory, no swap, a 50 GiB output allowance, and a 100 GiB free-disk gate.
The planning allowance is 0.5–8 hours on existing local resources. No GPU or
paid service is involved. Output is
`results/domains/full-candidate-architectures-v1/candidate_architectures.sqlite`.
The `queries` table stores candidate alternatives and policy uncertainty;
`proteins` preserves the taxon/protein-to-sequence mapping with foreign keys.

At launch, database construction is running. A successful producer receipt
will still require independent full readback against the source annotations,
competition outcomes and protein links. Domain architecture validation,
family integration, fusion/rearrangement inference and evolutionary tests
remain unfinished.

An independent full readback is queued under
`fungal-candidate-architecture-readback-20260922.service`. It waits for the
exact builder PID and creation time, requires a successful producer receipt,
and then streams every query against both the raw-hit database and the audited
competition table. `scripts/readback_candidate_domain_architectures.py` does
not import the architecture-producing function. It independently checks every
annotation field, coordinate ordering, model/type multiplicity, signature,
partial-HMM count, inclusive interval overlap, alternative reference and
policy-specific uncertainty count. It also compares all protein-link tuples
in sorted order and checks database integrity and foreign keys.

The validator passed 500 seeded four-policy cases and rejected six deliberately
corrupted exports (changed annotation type, repeat count, partial-match count,
uncertainty count, alternative reference and no-hit status). Reproduce with
`python scripts/check_candidate_architecture_readback.py`.
Its plan is `metadata/full_candidate_domain_architecture_readback_plan.json`;
state and the eventual completion receipt are under
`results/domains/full-candidate-architectures-readback-v1/`. Resources are one
CPU equivalent, 8 GiB memory, no swap and small audit outputs. The 0.5–8-hour
planning allowance starts after the builder completes. Queuing the audit does
not establish that the database has passed.

The independent audit has now passed in full. The archived receipt is
`metadata/full_candidate_domain_architecture_completed_readback.json`.
It checks 5,713,603 queries, 8,103,610 raw hits, and all 5,815,847 protein links
across 526 taxa. There are 53,248 queries whose retained hit sets differ among
the four competition policies. The 5,767,208 alternatives count sums alternatives
across queries; it is not a count of globally distinct biological architectures.
No-hit queries and uncertain overlap/nesting remain explicitly represented.

Family integration is the next stage; see
`docs/family-domain-bridge-20260922.md`. Passing this annotation audit does not
establish domain gains, losses, fusions, or rearrangements.
