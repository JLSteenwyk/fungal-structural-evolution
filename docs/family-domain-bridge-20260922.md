# Linking complete family partitions to candidate domain architectures

The bridge links every native gene ID to its taxon-specific protein ID,
sequence identifier and candidate-architecture source. It preserves both
expanded guide-specific family partitions: 658,183 profile-guide families and
658,522 MAFFT-guide families, each containing all 5,815,847 proteins from 526
taxa. These are family assignments, not completed reconciled orthology calls.

`scripts/build_family_domain_bridge.py` consumes the independently audited
architecture database and family partitions. It verifies the native ID files,
every family count and membership checksum, and the exact taxon/protein join.
Database uniqueness and foreign keys prohibit assigning a protein twice within
a guide or dropping a protein silently. The source family crosswalk remains
available in the `families` table; family names alone should not be compared
between guides.

The output is `results/domains/family-domain-bridge-v1/family_domain_bridge.sqlite`.
It contains `native_proteins`, `proteins`, `families`, and `assignments` tables.
To retrieve candidate alternatives, attach the pinned architecture database and
join its `queries` table to `proteins` using `sequence_id`. Annotation payloads
are not duplicated or reduced to a single preferred policy. The plan pins the
source database, receipts, native mapping, partitions, crosswalks, native MCL
reader and scripts. Inputs are rechecked before and after execution.

The pre-launch resource plan is `metadata/family_domain_bridge_plan.json`:
one CPU equivalent, 16 GiB RAM, no swap, 30 GiB output planning allowance and a
100 GiB free-disk gate. The 0.5–12 hour interval is uncalibrated planning, not
a measured ETA. The user service is `fungal-family-domain-bridge-20260922.service`.
It uses existing local CPU resources and does not restart GPU prediction.

Validation before launch:

```
.cache/envs/orthofinder/bin/python scripts/check_family_domain_bridge.py
```

The fixture verifies that identical protein identifiers in different taxa map
to their distinct sequences, both partitions retain both genes, and incomplete
domain mappings are rejected. A completed producer receipt still requires an
independent full bridge readback before use in evolutionary inference.
Domain-event inference, reconciliation and uncertainty-aware branch analyses
remain outstanding.

The initial independent full readback was queued under
`fungal-family-domain-readback-20260922.service`. It waits for the exact producer
PID and Linux process start ticks, then requires a successful producer receipt
and matching database checksum. `scripts/readback_family_domain_bridge.py`
does not import the producer or its native MCL reader: it reads the serialized
partition independently and checks every native ID, taxon/protein pair,
sequence/source link, family membership and provenance row. It also rejects
extra rows and checks SQLite integrity and foreign keys. Source and output
hashes are rechecked after the audit.

The readback plan is `metadata/family_domain_bridge_readback_plan.json`.
It requests one CPU equivalent, 8 GiB RAM, no swap, small receipt outputs and
an uncalibrated 0.5–12 hour runtime allowance after the producer exits.
`scripts/check_family_domain_bridge_readback.py` passed a complete fixture,
rejected changed sequence links and missing family members even after updating
the producer checksum, and rejected three malformed serialized partitions.
The initial output location was
`results/domains/family-domain-bridge-readback-v1/`; this attempt was later
stopped for the query-plan performance issue described below and did not pass.

Database construction has completed for both full partitions in 359 seconds.
The receipt is archived as `metadata/family_domain_bridge_completed_receipt.json`;
the independent full readback is now running. A downstream descriptive
within-family annotation inventory is queued behind that audit, documented in
`docs/family-architecture-variation-20260922.md`. No domain evolutionary event
is inferred from successful database construction.

## Indexed readback recovery

The first audit was live and consuming CPU, but SQLite chose the primary-key
index on `(guide, native_gene_id)` for each family-membership query, scanning
the guide's 5,815,847 assignments to filter one family. `EXPLAIN QUERY PLAN`
confirmed this access path. For the complete 16,019-member profile family
OG0000000, the original query took 12.26 seconds; explicitly selecting the
existing `(guide, family)` index took 0.0146 seconds and returned exactly the
same sorted membership hash. This is one query comparison under current cache
conditions, not a measured full-audit speedup or ETA.

The replacement `scripts/readback_family_domain_bridge_indexed.py` preserves
all checks and adds `INDEXED BY assignments_family` to that one query, plus
progress records at phase boundaries. It leaves the source database and its
indices unchanged. The full synthetic fixture passed; changed sequence links,
missing members and malformed partitions were still rejected. Reproduce that
check with `.cache/envs/orthofinder/bin/python
scripts/check_family_domain_bridge_indexed_readback.py`.

The old waiting architecture job was stopped first, then the old audit, with
process identities and state captured in
`metadata/family_domain_readback_performance_recovery.json`. Old scripts, plans
and incomplete outputs remain preserved. The replacement audit is running
under `fungal-family-domain-indexed-readback-20260922.service`, using
`metadata/family_domain_bridge_indexed_readback_plan.json` and fresh output
`results/domains/family-domain-bridge-readback-v2/`. It retains the one-CPU,
8-GiB, no-swap limits. The new downstream architecture job waits for this
exact process identity and a passing receipt. Full production validation is
still pending; the timing comparison is not its substitute.

The indexed production readback subsequently passed in 423.5 seconds. It checked
all 5,815,847 proteins and 526 taxa, all 658,183 profile-guide families and
658,522 MAFFT-guide families, complete source memberships and provenance,
SQLite integrity and foreign keys, and final input hashes. The receipt is
archived as `metadata/family_domain_bridge_completed_readback.json`, with its
producer, plan and script hashes checked on archival. The waiting architecture
analysis has started its profile-guide partition. This validates the bridge;
reconciled orthology and biological domain-event inference remain unfinished.
