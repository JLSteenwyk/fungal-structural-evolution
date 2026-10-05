# Full AFDB and ESMFold availability union

Every one of the **5,815,847 representative proteins across 526 taxa** has
been matched against the complete selected AFDB catalog and all completed
ESMFold models. All source model metadata, coordinate hashes, sequence
identities, lengths, provider/software and original prediction settings
remain distinct. No predictor is selected by comparing confidence scores.

The full independent reader reconstructs all **2,961,055 source models**,
every AFDB link and every representative FASTA sequence, emitted TSV and
SQLite row without importing producer logic. A separate FASTA parser is
used. All missing and overlapping source states, all 526 taxon rows and
the one ESMFold sequence outside the representative baseline are retained.
Both actual original tool waits, complete wrapper/native/manager transports
and **606 source/output bindings** close.
[Full completed union](../metadata/full_prediction_atlas_union_completed_20261005_v1.json)
and [independent source replay](../metadata/full_prediction_atlas_union_readback_20261005_v1.json).

| Existing model availability | Representative proteins |
| --- | ---: |
| AFDB only | 2,994,146 |
| ESMFold only | 24,801 |
| Both predictors | 722 |
| Neither source | 2,796,178 |
| At least one source | 3,019,669 |

Overall availability is **51.9%** before residue confidence/PAE qualification.
All **501 fungi and 25 outgroups** have at least one link. Fungi have
2,778,412 modeled representatives of 5,431,687 (**51.2%**); outgroups have
241,257 of 384,160 (**62.8%**). ESMFold adds 24,163 fungal and 638 outgroup
links beyond AFDB. Both-source coverage occurs only in fungi in this frozen
baseline. Taxon representation does not establish equivalent protein
coverage across lineages or eliminate ascertainment bias.
[Complete 526-row source table](../metadata/full_prediction_atlas_union_taxon_coverage_20261005_v1.tsv).

This step integrates **existing** predictions; it inferred no new structures.
The 2,796,178 representatives without a model in these sources remain in
scope. They are not assumed to lack models elsewhere, lack homologs or
represent losses or novelty. The domain registry currently joins the
[AFDB source](completed-retrieval-domain-registry-20261005.md); ESMFold
domain integration and full residue-confidence/PAE auditing remain required.

## Original source custody and execution

The union retains **2,935,733 AFDB models and 25,322 ESMFold models**. All
seven ESMFold prediction configurations remain explicit. Original local
representation versions are not AFDB release identifiers. Full raw model
metadata are retained in `model_metadata.jsonl`, including original PAE
URLs/flags or local NPZ paths; a snapshot flag is not a fresh filesystem
availability audit. All 25,322 ESMFold coordinate files, PAE NPZ files and
prediction receipts were rehashed: **7,339,969,721 coordinate bytes** and
**25,846,628,960 PAE bytes**. AFDB coordinate hashing relies on the recent
complete 1.10 TB catalog audit and its bound model table; it is not repeated
by the union. The independent reader does not repeat those byte audits.

Producer original tool session **58078** returned zero after 8 minutes
41.507 seconds. The independent reader's original session **4670** returned
zero after 7 minutes 21.008 seconds. Both were launched only after required
predecessor closure. The producer used two CPUs, 32 GiB RAM, no swap,
28 GiB address space and a 64 GiB output allowance. The reader used two CPUs,
64 GiB RAM, no swap and 56 GiB address space. Both used one BLAS thread;
14,400 CPU-second and 21,600 wall-second caps are not ETAs. Their prelaunch
0.25–3 hour ranges were uncalibrated. No GPUs, new predictions, downloads,
paid infrastructure or charges were used.
[Producer resources](../metadata/full_prediction_atlas_union_resources_20261005_v1.json)
and [reader resources](../metadata/full_prediction_atlas_union_readback_resources_20261005_v1.json).

## Annotated product outside the representative baseline

One ESMFold sequence has no exact representative link. It is the 769-residue
annotated protein `KAF8745614.1` in *Amanita brunnescens* (`F87326`), marker
`5006412at2759`. The frozen mapping assigns it to `gene-AX14_006929`; the
baseline selected the longer 1,072-residue annotated product `KAF8745613.1`
from the same gene. The longer product has neither source model in this
union. The shorter model remains in the model inventory with its original
configuration and source annotation. The baseline has not been changed to
obtain a match. Expression, biologically correct isoform assignment and
gene duplication are not established by these annotations.
[Complete exception diagnostic](../metadata/full_atlas_unlinked_esmfold_20261005_v2.json).

The first diagnostic asserted that a single-marker FASTA hash equaled the
full QC proteome hash; those are different files. The failed code, original
exit-one tool result, complete wrapper logs and failure audit are preserved.
The new diagnostic binds both files separately and compares the marker
sequence directly to its full source product and ESMFold hash, and the longer
product directly to the baseline hash. The biological sequences match.
[Retained first failure](../metadata/full_atlas_unlinked_esmfold_failure_audit_20261005_v1.json).
The diagnostic's provisional readback flag describes its own limited scope;
the authoritative full union/readback completion is the final receipt above.

## Reproduction and remaining work

Large outputs remain outside Git at
`results/structures/full-prediction-atlas-union-20261005-v1/`:

- `prediction_atlas.sqlite` — all model, original AFDB link and representative records
- `protein_source_links.tsv` — every representative, including neither-source rows
- `model_metadata.jsonl` — complete original model records, with source retained
- `taxon_source_coverage.tsv` — every taxon denominator and source disposition

[The pinned plan](../metadata/full_prediction_atlas_union_plan_20261005_v1.json)
binds the complete catalog, representative inventory, sampling manifest,
original ESMFold source/readback/configuration chain and code. Use fresh
output/receipt paths and rebuilt input pins for an independent reproduction.
`build_full_prediction_atlas_union_v1.py` produces the complete union;
`readback_full_prediction_atlas_union_v1.py` requires its actual original
execution transport and reconstructs all sources/records;
`close_full_prediction_atlas_union_v1.py` binds completion and releases a
byte-identical taxon table. Cloning the repository does not supply the
large raw inputs; their public release remains outstanding.

Next audit confidence and PAE on the complete source/model identities,
integrate ESMFold into domain/structural-family workflows, expand source
controls without treating raw structural distances as additive, and retain
all missingness and isoform exceptions in phylogenetic analyses. Accepted
phylogenetic/reconciliation/dating frameworks, adequate ancestral uncertainty,
calibrated effects and all eight evolutionary aims remain incomplete.
