# Whole-proteome structure catalog

The marker catalog does not measure whole-proteome coverage. A separate
full-data catalog has completed screening all **5,815,847 representative proteins
from 526 sampled entries**. It uses the same frozen AlphaFold retrieval log
as the refreshed marker work and preserves the one-longest-protein-per-gene
selection policy of `metadata/gene_representatives_receipt.json`.

`scripts/catalog_whole_proteome_structures.py` replays the last recorded
status of each UniProt accession, retains verified GDM/AlphaFold Monomer
v2.0 pipeline models, and selects one model per exact sequence hash by
mean C-alpha pLDDT, model version and model identifier. It independently
reads and hash-checks every representative FASTA, checks sequence length,
preserves every species/protein link even when sequences are identical,
and hash-checks every selected coordinate file against its original
verification record. It does not reparse CIF content or claim confidence
qualification. Source selection is separate from the empirical comparison
of prediction methods required downstream.

The immutable output directory is
`results/structures/whole-proteome-afdb-catalog-20260922-v1`. Outputs include
all protein/model links, model provenance in JSONL, per-taxon coverage and
a completion receipt. Downstream work must require that completion receipt;
partial files alone do not establish successful completion. Absence from this catalog
means absence from the frozen source snapshot, not absence from public
databases or from later downloads. ESMFold models are not included in this
source-specific catalog and will require an explicit subsequent integration.

The plan is `metadata/whole_proteome_structure_catalog_plan.json`, launched
as `fungal-whole-proteome-structure-catalog-20260922.service`. Resources are
one CPU equivalent, 32 GiB RAM, zero swap, reduced CPU/IO scheduling priority,
5 GiB output allowance and an uncalibrated 0.5–12 hour runtime range. The
filesystem had approximately 12 TiB free before launch. All coordinate
access is read-only; no GPU prediction or paid resources are used.

`python scripts/check_whole_proteome_catalog.py` exercises full CLI behavior
on a synthetic fixture: best-confidence model selection, latest failed
accession exclusion, wrong-provider exclusion, repeated sequence links
and rejection of a changed coordinate file. Fixture coordinate bytes are
not valid structures: this explicitly tests catalog and integrity logic,
not coordinate validation. The scientific catalog and independent full-data link readback have passed.

Reproduce with a fresh output path in a new pinned plan:

```bash
python scripts/catalog_whole_proteome_structures.py \
  --plan metadata/whole_proteome_structure_catalog_plan.json
```

Confidence/PAE qualification, domain-level comparisons, structural clusters,
remote-homology review and atlas-wide evolutionary tests remain unfinished.

## Independent full-data link readback

`scripts/readback_whole_proteome_catalog.py` completed after waiting for the exact
producer PID and creation time. It requires successful catalog completion,
checks the frozen primary inputs and output hashes, then independently
replays latest accession statuses and selects models by a stable descending
sort. A separate FASTA parser reconstructs every expected link across all
5,815,847 proteins and compares it against the emitted table. It verifies
all per-taxon coverage rows, model provenance, totals and absence of extra
links. It does not repeat the producer's coordinate-byte hashing or perform
new CIF-content/PAE validation.

The readback plan is
`metadata/whole_proteome_structure_catalog_readback_plan.json`, launched as
`fungal-whole-proteome-catalog-readback-20260922.service`, with one CPU,
32 GiB RAM, no swap and reduced scheduling priority. It writes only a small
receipt, with an uncalibrated 0.5–6 hour planning range after the producer
finishes. The updated synthetic fixture verifies both successful readback
and rejection of an incorrect protein link whose artifact checksum has
been updated. This tests semantic reconstruction rather than only hashes.

The completed catalog contains **1,290,278 unique models linked to 1,319,513
of 5,815,847 proteins (22.69%)**, with models for 496 of 526 sampled entries.
The producer verified 492,263,006,636 coordinate bytes. Independent readback
reconstructed all representative sequences, model selections and protein links.
Both systemd jobs terminated successfully. Completion receipts and the complete
526-row coverage table are archived under `metadata/whole_proteome_structure_catalog_completed_receipt.json`,
`metadata/whole_proteome_structure_catalog_completed_readback.json` and
`metadata/whole_proteome_structure_taxon_coverage.tsv`. All catalog artifact
hashes and coverage-table totals were rechecked before archival.

| Sampling role | Entries with models / total | Proteins with models / total | Coverage |
|---|---:|---:|---:|
| Fungal ingroup | 472 / 501 | 1,254,497 / 5,431,687 | 23.10% |
| Non-fungal outgroup | 24 / 25 | 65,016 / 384,160 | 16.92% |

These are sampled entries, not a new validation of distinct species identities.
Coverage is source-specific availability before confidence filtering. It does
not include local ESMFold models and does not establish that proteins missing
from this catalog lack models elsewhere.

## Search-database construction running

The next stage passed the independent-readback gate and is running as
`fungal-whole-proteome-foldseek-db-20260922-v2.service`. Its plan,
`metadata/whole_proteome_foldseek_database_plan.json`, pins the catalog and
readback plans, installed executable and implementation. It will process
the entire verified catalog into
`results/structural_clusters/whole-proteome-afdb-database-20260922-v1`.

The installed Foldseek version e3fadcd07f971e864c094ac4f3a78bf4ed845e07
accepts a newline-separated path list, verified against its [pinned source](https://github.com/steineggerlab/foldseek/blob/e3fadcd07f971e864c094ac4f3a78bf4ed845e07/src/strucclustutils/structcreatedb.cpp)
and a native two-model compatibility fixture. This avoids creating a million
additional symlinks. The command explicitly sets four threads, GPU off,
pLDDT70 seeding masking and float32 coordinate storage.

After native database construction, every lookup identity and amino-acid
sequence hash is checked against the catalog. The verifier checks every
3Di record's length and alphabet and every coordinate record's size and
finiteness. It does not independently regenerate 3Di values or compare
all coordinates back to CIF. Corrupted amino-acid and nonfinite-coordinate
fixtures were rejected. pLDDT seeding masks are not full confidence masks.

Resources are four CPU equivalents, 64 GiB RAM, no swap, 250 GiB output
allowance, 500 GiB free disk and 128 GiB available-memory gates. The
uncalibrated runtime allowance is 1–48 hours after catalog verification.
No pairwise search or clustering is launched by this database stage;
those require their own measured dimensions, resource estimates and
membership/edge checks. The launch receipt also records an initial setup
failure before database work began, followed by the verified corrected
launch. No existing database or running analysis was restarted.

## Gene-family coverage bridge completed; independent readback pending

`scripts/bridge_whole_proteome_structures_to_families.py` completed after
the independent catalog readback passed. Its plan is
`metadata/whole_proteome_family_coverage_plan.json`, launched as
`fungal-whole-proteome-family-coverage-20260922.service`. The stage requires
the previously completed independent family/domain bridge audit, then joins
every structural link by exact taxon, protein accession and sequence hash.
The existing family database is attached read-only; results are written to
`results/structures/whole-proteome-family-coverage-20260922-v1`.

Both complete partitions are retained: 658,183 profile-guide families and
658,522 MAFFT-guide families, each covering 5,815,847 proteins. Output rows
include total family proteins and taxa, proteins and taxa with models,
distinct modeled sequences and distinct models. Families with no structures
remain present. Same-taxon copies and identical sequences are not confused
with independently sampled taxa or models. A hand-calculated fixture checks
these distinctions and the retention of an entirely unmodeled family.

The stage uses one CPU equivalent, 32 GiB memory, no swap, 20 GiB output
allowance and a 100 GiB free-disk gate, with an uncalibrated 0.5–24 hour
runtime range. Both partition totals and every structural link's sequence
identity must agree before completion. The new coverage outputs will still
require independent readback. These family assignments are not finalized
reconciled orthology, and the models have not been confidence-qualified
atlas-wide. Coverage does not establish duplication or domain events,
remote homology, structural acceleration or statistical power.

The producer reports 30,959 profile-guide families and 30,934 MAFFT-guide
families with models in at least two taxa. All 1,319,513 modeled protein links
occur in each complete family partition. These are provisional producer totals:
independent family-level reconstruction is the next acceptance gate. A family
with every protein modeled may be a singleton and is not necessarily eligible
for comparative analysis. No family-level evolutionary result is inferred here.

## Independent family coverage verification running

`scripts/readback_whole_proteome_family_coverage.py` now checks both complete
family partitions under `fungal-whole-proteome-family-readback-20260922.service`.
Its plan, `metadata/whole_proteome_family_coverage_readback_plan.json`, pins
15 source, implementation and output files. The checker reconstructs the
native gene/model identities from the catalog and original protein database,
then compares every identity against the new structural bridge. It separately
streams each source family and accumulates Python sets for taxa, sequences
and models; it does not reuse the producer's SQL coverage aggregation.
Every coverage-table row, including families without models, must match.

The checker also reports counts of families with models in at least 2, 4,
10, 25, 50 and 100 taxa. These describe availability, not reconciled orthology,
independent evolutionary transitions, confidence qualification or statistical
power. The hand-calculated fixture includes same-taxon copies, shared model
sequences and an unmodeled family; incorrect model identity, sequence identity,
taxon counts, extra families and omitted unmodeled families are rejected.
Reproduce those checks with
`python scripts/check_whole_proteome_family_readback.py`.

The full run uses one CPU equivalent, 32 GiB RAM, no swap, reduced scheduling
priority and a 0.01 GiB output allowance. The pre-launch runtime planning range
is 0.1–12 hours and is uncalibrated. The small final receipt will be written to
`results/structures/whole-proteome-family-coverage-20260922-v1-readback.json`;
launch and fixture records are versioned in `metadata/`. The scientific
family-level verification remains pending until that receipt passes.
