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
