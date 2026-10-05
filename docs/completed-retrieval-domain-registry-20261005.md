# Full refreshed structural-domain registry

The [complete catalog and independent change replay](completed-retrieval-catalog-refresh-20261005.md)
now feed a fresh annotation registry across **all 2,935,733 selected AFDB
models and 2,994,868 representative-protein links**. This advances the domain
atlas used by domain evolution and duplication analyses. Full producer,
independent reconstruction and actual original execution closure now pass.
[Completed registry](../metadata/completed_afdb_domain_registry_completed_20261005_v1.json).
The fixed study
universe remains 501 fungi and 25 outgroups; source-unlinked proteins stay
explicit in the full coverage table rather than being interpreted as losses.

The qualified `build_structure_domain_registry.py` is unchanged. It joins
each exact sequence to the existing complete annotation database, retains
alignment and envelope positions, all Pfam types and all four architecture
policies. Candidate Domain intervals retain the original conservative
architecture, HMM coverage of at least 0.70 and minimum alignment length of
30 residues. This is an annotation interval registry, before residue
confidence/PAE, validated structural boundaries, homology or evolutionary
event inference. ESMFold remains a separate prediction source for subsequent
integration; no model or architecture rule is changed.

[Immutable full producer plan](../metadata/completed_afdb_domain_registry_plan_20261005_v1.json)
and [prelaunch resource estimate](../metadata/completed_afdb_domain_registry_resources_20261005_v1.json)
bind the full closed catalog, complete 11.72 GB annotation database and
original code. The new adapter requires both catalog and comparison closure
before invoking the unchanged producer. A separate reader keeps the complete
qualified reconstruction of every model, interval, policy membership,
candidate flag and protein link. Its new gate requires the actual original
producer wait/transport, complete source pins and receipt hash; PID
disappearance alone is insufficient.
[Reader source equivalence](../metadata/completed_afdb_domain_registry_source_equivalence_20261005_v1.json)
and [full reader plan](../metadata/completed_afdb_domain_registry_readback_plan_20261005_v1.json).

The original producer was wrapper PID 251281, created at Unix time
1791219275.13, invocation
`544cc1e796ce41658556ca4aca061ef7`, original tool session **93892**.
Its actual original tool wait returned zero after 9 minutes 5.495 seconds.
[Original launch identity](../metadata/completed_afdb_domain_registry_launch_20261005_v1.json).
The independent reader's original session **49632** returned zero after
10 minutes 21.618 seconds. It reconstructs every model, interval, policy,
candidate flag and protein link. Both original whole journals and all
56 source/output bindings close. The separate final all-taxon aggregation
and closure ran in 49.907 seconds under two CPUs and 16 GiB/no swap.

The registry contains **2,944,169 unique annotation intervals** and
**11,718,744 policy memberships**. Counts of candidate Domain intervals
depend on the original policy:

| Architecture policy | Candidate intervals per model |
| --- | ---: |
| Alignment bitscore | 1,351,205 |
| Alignment E-value | 1,351,926 |
| Envelope bitscore | 1,347,031 |
| Envelope E-value | 1,347,757 |

Policies are alternatives; counts must not be added as independent evidence.
The [complete 526-row taxon table](../metadata/completed_afdb_domain_registry_taxon_coverage_20261005_v1.tsv)
retains all 5,815,847 denominators and 2,820,979 AFDB-unlinked proteins. It
counts interval occurrences per linked protein, which can exceed model
interval counts when an identical model is reused. Under alignment E-value,
1,000,200 linked proteins have at least one candidate interval, with
1,384,603 total candidate occurrences. There are 1,930,398 linked proteins
with raw Pfam hits and 1,064,470 without; absence of Pfam hits is not lack
of function or novelty. No physical boundaries or evolutionary events are
accepted by these annotations.

The producer uses two CPU equivalents, 32 GiB RAM, no swap, 28 GiB address
space, one BLAS thread, a 32 GiB output allowance, 7,200 CPU-second and
14,400 wall-second caps. Its 0.1–1.5 hour planning range is uncalibrated.
The previous 1,910,138-model registry took 347.791 seconds and generated
a 2.018 GiB SQLite database. The new scope is 1.537 times that model count;
disk contention and annotation complexity limit timing extrapolation.
The full reader has separate two-CPU/16 GiB/no-swap caps and a 12 GiB address
space limit. No annotation search, prediction, GPU use or new cost is launched.

Outputs remain outside Git at
`results/domains/whole-proteome-structure-domain-registry-20261005-v1/`
and the completed `...-readback/` directory. Both final receipts, full source
reconstruction and actual original execution closure are now available.
Do not restart or overwrite these immutable outputs. For an
independent reproduction, create new output/receipt paths and rebuild the
complete input pins before running the adapter or reader.

The full structural atlas, accepted phylogenetic framework, reconciliations,
adequate ancestral uncertainty and all eight evolutionary aims remain
required and incomplete.
