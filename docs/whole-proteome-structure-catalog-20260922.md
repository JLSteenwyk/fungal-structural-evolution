# Whole-proteome structure catalog

The marker catalog does not measure whole-proteome coverage. A separate
full-data catalog is running against all **5,815,847 representative proteins
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
files from the running producer are incomplete. Absence from this catalog
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
not coordinate validation. The scientific catalog remains in progress.

Reproduce with a fresh output path in a new pinned plan:

```bash
python scripts/catalog_whole_proteome_structures.py \
  --plan metadata/whole_proteome_structure_catalog_plan.json
```

Confidence/PAE qualification, domain-level comparisons, structural clusters,
remote-homology review and atlas-wide evolutionary tests remain unfinished.

## Independent full-data link readback

`scripts/readback_whole_proteome_catalog.py` is queued behind the exact
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

The producer's live log reports 1,290,278 selected models and has entered
coordinate hash verification. This is a provisional selection count,
not a completed coverage result; acceptance awaits both completion receipts.

## Search-database construction queued

The next stage is queued behind the independent readback as
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
