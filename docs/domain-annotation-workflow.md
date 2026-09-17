# Domain annotation workflow

The first production layer annotates all 59,840 complete marker-protein records spanning the 526-taxon extraction. Exact sequence deduplication yields 58,883 search targets (33,047,630 residues), with all marker, taxon, protein and source checksums retained in a link table. Deduplication reduces repeated search work; it does not collapse species-specific records for evolutionary analysis. Full representative-proteome annotation remains required and is a subsequent stage.

## Pinned resources and execution

Pfam 38.2 contains 30,134 families and is based on UniProtKB 2025_03 according to its version file. The [release-specific EBI archive](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/) supplies profile HMMs, entry types, clan assignments and active-site metadata. `prepare_pfam.py` pins publisher MD5 values, verifies decompression and records local SHA256 hashes. Active-site metadata are acquired for later use; no site transfer is implemented or claimed here. All HMMs must have unique accessions, match-state lengths and gathering thresholds. Their total count must agree with the release. Profiles are assigned deterministically to 64 chunks balanced by summed match-state lengths.

The environment is specified in `environments/domains.yml`; the current installed HMMER binary is additionally pinned by hash at execution. HMMER 3.4 `hmmsearch --cut_ga --cpu 2 --noali` searches each HMM chunk against the entire deduplicated marker database. Four jobs run concurrently. E-values refer to this unique-sequence target database, not all fungal proteins; Pfam-specific sequence and domain gathering scores control acceptance. See the [Pfam glossary](https://pfam-docs.readthedocs.io/en/latest/glossary.html) for the two curated gathering thresholds. Search commands, binary hash/version, source hashes, runtimes and output hashes are recorded. Resume accepts only matching completed receipts; unreceipted outputs require inspection after interruption. A lock prevents concurrent writers.

```bash
python scripts/prepare_pfam.py
python scripts/prepare_marker_domain_inputs.py
python scripts/run_marker_domains.py --output results/domains/marker-search-v1
python scripts/summarize_marker_domains.py --search results/domains/marker-search-v1 --output results/domains/marker-annotations-v1
```

Marker inputs and interpreted annotations are immutable snapshots; prepare them in a new versioned location when inputs change. The search can resume its existing output directory under an unchanged configuration. Raw files remain under ignored `data/pfam/`, `data/domains/` and `results/domains/`. Published metadata receipts contain retrieval and verification information.

## Interpretation and downstream gates

Only a complete search across all pinned profiles can enter the annotation summary. The parser verifies target lengths, model identities and lengths, inclusive alignment/envelope coordinates, finite scores and gathering-threshold consistency (allowing 0.051 bits for printed-score rounding). It preserves Pfam versioned accession, entry type, clan, HMM coverage, scores, E-values and coordinates. Types such as Family, Repeat and Domain remain distinct; a Pfam match is not automatically a discrete structural domain or an orthogroup.

All gathering-threshold matches are retained. Overlaps in alignment coordinates, including a shared boundary residue, are explicitly recorded with same-clan and same-family flags. This is not a resolved domain architecture. Overlap resolution, fragmented matches, nested domains, repeated domains, alignment-boundary uncertainty and independent structural support must be assessed before counting gains, losses, fusions or duplications. Unknown or undetected domains are not treated as confirmed absences. Preserving the raw hits supports later comparisons of alternative overlap policies.

The parser and overlap-coordinate tests are implementation checks, not biological validation. Domain-restricted structural comparisons, domain architecture reconstruction, functional-site transfer and full-proteome annotation remain pending until their execution is separately documented.

## Full representative-proteome inputs

`prepare_full_domain_inputs.py --output data/domains/full-inputs-v1` prepares the complete 526-taxon representative baseline for domain annotation. Each source FASTA is checksum-verified against the completed representative receipt; unique protein IDs and per-taxon counts must match. Exact sequence matches to the marker-domain input retain links to that query set, and each additional exact unique sequence is written once. `protein_links.tsv` retains every selected taxon/protein ID and identifies the query source. This preserves species-specific and gene-specific associations even when the search is deduplicated. Unresolved gene mappings remain unresolved and cannot become inferred gene-copy counts through this preparation.

The marker and additional search partitions have different target database sizes. Curated gathering scores can be applied consistently, but raw sequence/domain E-values must remain labeled by search partition; conditional domain E-values cannot be globally calibrated by naively multiplying them by a protein-count ratio. Combining hits requires verified completion against the same Pfam version and compatible HMMER settings. Input preparation alone does not complete full-proteome annotation, and marker results cannot be reused before their complete search passes validation.

The immutable input receipt records source/manifest hashes, per-taxon counts, unique sequence/residue totals and output checksums. Large FASTA and link tables remain outside Git. Its resource estimate precedes launch in `docs/resources.md`; the full-proteome search requires a separate estimate based on observed production throughput.

## Queued full-proteome execution

The additional full-proteome search is queued to reuse the marker search's CPU allocation:

```bash
flock results/domains/marker-search-v1/.lock bash -c 'python scripts/summarize_marker_domains.py --search results/domains/marker-search-v1 --output results/domains/marker-annotations-v1 && python scripts/run_full_domains.py --output results/domains/full-search-v1'
```

This command waits for the live marker process to release its lock, then requires a complete, validated marker annotation before launching four two-worker-thread HMMER jobs against the additional sequences. The full runner independently checks exact marker profile coverage, raw table and annotation hashes, marker/full-input linkage, identical Pfam receipt and identical HMMER binary. Failed or incomplete marker results stop the dependency chain. The new runner preserves per-chunk verified resume and streams checksums; output-table completion checks read only the file tail to avoid loading large result files into memory.

`metadata/full_domain_execution_plan.json` pins the queued command, scripts, input receipt and resource estimate at submission. It is a historical queue record, not a live status or completion receipt. The queue holds the marker lock through the dependent full search; do not launch a second competing queue or restart either search merely because a polling interval expires. Once marker annotations exist, a future resume should invoke `run_full_domains.py` directly under its own output lock; do not rerun the immutable annotation writer. Merging the two search partitions and resolving architectures remain subsequent steps.

## Tree mapping of detectable Pfam-profile changes

`results/domains/marker-pfam-detection-changes-v1` maps detectable profile states to the 70 audited marker trees in snapshot v5. Domain, Family and Repeat types remain separately labeled. All 1,310 marker/profile combinations are retained under two policies (2,620 disposition rows). For the baseline, any gathering-threshold hit is detected; no hit is undetected. The conservative policy treats a protein with any overlapping annotated hits as unknown for every profile and treats a profile with only hits below 70% HMM coverage as unknown. This deliberately conservative masking is not an overlap-resolution algorithm. Neither policy converts no detectable hit into confirmed biological domain absence.

At least five detected and five undetected sampled tips are required for mapping. This leaves 185 combinations in 57 markers under the baseline and 81 in 44 markers under the conservative policy. Full source sequence identities are retained, tree-tip membership is exact, and observed residues in each trimmed marker alignment are checked for preserved order in the corresponding domain-query sequence. Profiles failing eligibility remain visible rather than being silently omitted.

A binary symmetric unit-cost Sankoff calculation finds the minimum number of detection-state changes on each unrooted marker tree. Outside-subtree costs identify edges where a change occurs in at least one optimal reconstruction or in every optimal reconstruction. No single ancestral reconstruction is selected. The 6,631 possible-change edge rows retain marker splits, terminal/internal status and reported SH-aLRT support. These are conditional minimum-change locations, not calibrated probabilities, supported species-lineage events or rooted gains/losses. Missing-marker bias, gene-tree uncertainty, undetected profiles, sequence fragments, nested domains and competing annotations remain material limitations. Copy counts, order changes, fusions, rearrangements and reconciled gene-family turnover still require separate analyses.

The algorithm matched exhaustive enumeration for all 162 leaf-state cases on two small tree representations, including unknown states, checking scores and every possible/mandatory edge. An independent empirical readback rebuilt all 2,620 dispositions and 126,667 eligible tip states from raw annotations and recomputed all 266 scores with set-based Fitch traversal. Empirical edge marginals were checked for counts/uniqueness and the mandatory-edge bound, not independently exhaustively enumerated. Among the 81 combinations eligible under both policies, 46 minimum scores decrease after ambiguity masking; none increases, as required when allowed tip states are enlarged. Median scores of nine versus eight across the different eligible sets are descriptive and should not be interpreted as a paired effect.

Scripts are `map_marker_pfam_detection_changes.py` and `readback_marker_pfam_detection.py`. Versioned receipts, profile dispositions, resource estimates and verification records are under `metadata/marker_pfam_detection_*`. Full tip and possible-change tables stay outside Git. Planning used one CPU, 4 GiB memory, 2 GiB output and 0.02–2 hours on the existing host. Full-proteome Pfam searches continue; this marker-based analysis does not substitute for full-family architecture reconstruction.

## Streaming annotation of the additional full-proteome search

The full-search producer has 52 completed raw chunks and 12 pending out of 64. An immutable snapshot of those 52 receipts (6,498,460 raw hits) is pinned in `metadata/full_pfam_annotation_snapshot_plan.json`; actively written tables are excluded. `annotate_completed_full_pfam_chunks.py` validates all input artifact hashes and the full additional query FASTA (5,654,720 sequence hashes; 2,520,735,813 residues), then streams each completed raw table through the previously validated Pfam parser. Query identity, profile identity, sequence/HMM coordinates, finite values, gathering thresholds and row counts must match. Each output shard has its own source/output hashes and stable chunk/row hit ID.

The outputs preserve Domain/Family/Repeat and other Pfam types, clan/description metadata, partial hits, scores, conditional and independent E-values. Every row is labeled `additional_full_proteome`; E-values are not rescaled or pooled with the smaller marker search. The snapshot remains incomplete until all 64 search chunks are represented. It does not resolve overlaps or establish absence for unannotated proteins.

The live transient service is `fungal-full-pfam-annotation-20260916.service`, initially PID 321008; verify its recorded creation time rather than assuming that PID remains valid. It uses one CPU, an 8 GiB planning allowance and a 16 GiB service limit, with 10 GiB output allowance and 0.1–6 hours planned. Available-memory and disk gates passed before launch. Source scripts are pinned and must remain unchanged while running. At the recorded checkpoint, 18 shards with 2,472,034 rows had completed; the full snapshot was not yet complete.

A separate live waiter (identity in `metadata/full_pfam_annotation_readback_launch.json`) runs `readback_full_pfam_annotation_chunks.py` only after this producer exits and the completed annotation receipt is available. It independently compares every annotated row with the original raw field order: identifiers, coordinates, scores, E-values, posterior accuracy and HMM coverage. Descriptive Pfam metadata is supplied by the pinned parser and is not independently parsed in this readback. Planning is one CPU, 2 GiB memory and 0.05–3 hours. Failed or partial snapshots cannot pass as complete. The waiter and service do not survive a host reboot.

After successful readback, merging the remaining chunks and the reusable marker partition, retaining every representative-protein link, detecting overlaps and reconstructing architecture uncertainty remain required. Annotation of the completed shards is production progress across the full proteome dataset, not a substitute for those steps.

### Remaining-shard completion handoff

The initial 52-shard annotation snapshot completed all 6,498,460 rows; its independent raw-field readback is in progress. `fungal-full-pfam-completion-20260916.service` waits for both the original full-search producer and this readback producer using PID, creation time and exact command identity. Successful prerequisite receipts are required after they exit. It then creates a source-pinned plan for exactly the remaining 12 chunks, annotates them, performs the same full raw-field readback, and combines the two disjoint audited snapshots into `full-annotation-catalog-v1`.

The catalog gate requires all expected Pfam chunks exactly once, matching source configuration and profile hashes, matching shard checksums and raw hit counts, and full profile/hit totals. Synthetic checks verified acceptance of a complete disjoint fixture and rejection of a missing chunk, duplicate chunk and changed annotation artifact; these are orchestration checks, not empirical validation. Source snapshot files remain separate and unchanged; the catalog records their paths and hashes without duplicating large tables.

The service is pinned by `metadata/full_pfam_annotation_completion_config.json`; its script is `advance_full_pfam_annotation.py`, with catalog assembly in `combine_full_pfam_annotation_snapshots.py`. It uses one CPU, an 8 GiB planning allowance, 16 GiB service limit, 32 GiB available-memory/50 GiB free-disk gates and a 0.1–9 hour planning range after search completion. No new charges. Failure leaves reviewable state and does not automatically restart or reinterpret incomplete output as success. A transient service does not survive reboot. The marker-query partition, protein-level joins, overlap resolution and architecture inference remain separate downstream requirements.

The initial raw-field readback subsequently passed for all 52 chunks and all 6,498,460 rows. Its receipt is archived as `metadata/full_pfam_annotation_readback_receipt.json`. The completion controller continues waiting for the remaining HMMER search; no full64 catalog exists yet.

## Full-proteome query linkage readback

`scripts/readback_full_domain_protein_links.py` independently reconstructs every query identifier from source representative FASTA sequences, using a minimal streaming parser separate from the preparation script's Bio.SeqIO parser. It requires an exact match for every taxon, protein accession, sequence hash and search partition in `protein_links.tsv`. It checks unique protein identifiers within each taxon, all query sequence hashes, disjoint marker/additional query sets, complete additional-query coverage, and per-taxon and aggregate counts against pinned source receipts. Marker queries not represented in the full proteome are reported explicitly. Identical sequences remain linked to every source protein; deduplication is only a search optimization.

```bash
python scripts/readback_full_domain_protein_links.py \
  --output results/domains/full-protein-link-readback-v1
```

The output is immutable; use a fresh directory for subsequent runs. A passing receipt verifies computational linkage, not the biological correctness of representative selection, taxonomy, homology or Pfam assignments. Downstream joins must retain the query-source partition and pass both this linkage check and the corresponding complete annotation checks before interpreting protein-level coverage. Full architecture resolution and evolutionary inference remain separate requirements. The resource estimate is recorded in `metadata/full_domain_protein_link_readback_plan.json`.

The full readback completed on September 16: all 5,815,847 protein links across 526 taxa passed. It independently verified 5,713,599 distinct full-proteome sequences, including 58,879 reused marker queries with 60,033 protein links; four marker queries are outside the representative-proteome universe. All 5,654,720 additional queries are represented, with no overlap between query partitions. The receipt and 526-row taxon summary are archived under `metadata/full_domain_protein_link_readback_*`. Synthetic re-pinned fixtures confirmed rejection of missing links, extra duplicate links, wrong partitions and wrong taxa. This establishes linkage integrity; complete annotation merging and architecture analysis remain pending.


## Complete partition and protein joins (2026-09-17)

The full raw-hit database is being built by
`scripts/build_full_domain_database.py` under
`fungal-full-domain-database-20260917.service`. All 7,939,960 additional hits
and 163,650 marker hits are retained alongside the 5,815,847 representative
protein links. The query universe includes 5,713,599 full-proteome sequences
plus four marker-only sequences; those four have no representative-protein
link. Exact source hit IDs remain qualified by search partition. Coordinates,
scores and partition-specific E-values are stored as their original text;
no numerical normalization, overlap filtering or architecture assignment is
performed. Consumers must explicitly cast numerical fields for numeric SQL
comparisons rather than use text sorting.

The database exposes `protein_hits` for hit joins and `protein_detection` for
all proteins, including the explicitly qualified `no_GA_hit_not_proven_absence`
state. Identical sequences share search results while every species/protein
association remains available. Unique protein keys, unique partition/hit IDs,
query partition integrity and complete counts are checked during construction.

After construction, `scripts/readback_full_domain_database.py` compares every
stored hit field and every protein link to the immutable source tables. A
separate source-hit counter supplies all 526 expected taxon-level detection
and joined-hit totals, which are compared with database aggregation. Its
receipt and `taxon_domain_detection.tsv` are required before using the joins.
No database completion or readback is claimed yet. Clan-aware overlap review,
nested domains, architecture uncertainty and family-tree integration remain
subsequent analyses.

The plan and exact live launch identity are versioned as
`metadata/full_domain_database_{plan,launch}.json`. The existing host supplies
one CPU, an 8 GiB memory planning allowance and 16 GiB service cap; 100 GiB
free disk is required and the build has a 50 GiB database allowance. The
0.25–6 hour window is a broad planning estimate. This transient service does
not survive reboot; partial outputs require review and are not overwritten.


## Preserve curated nested-domain relationships

Pfam describes nested domains as insertions within a containing domain and
marks them with NE tags ([Pfam database paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC308855/)).
The checksum-verified 38.2 metadata contains 160 NE occurrences representing
150 distinct directed model relationships in 111 containing families. Twenty-seven
families have multiple NE tags (maximum11), and ten occurrences repeat a
pair already recorded in that source. `extract_pfam_nested_domains.py` retains
every occurrence and resolves model names to accessions within the same release.
An independent regular-expression/block parser confirmed the complete exported
pair multiset, including multiplicity. Receipts and the small relationship table
are versioned as `metadata/pfam_nested_*`.

The existing marker annotation parser retains one value per metadata key; it
is sufficient for its scalar hit annotation fields but must not be reused to
extract this multivalued NE relation. Its existing annotations are unchanged.
Architecture processing must use the complete nested relationship table and
observed coordinates; these curated model relationships alone do not establish
nesting in any project protein or justify accepting all overlapping hits.


The full database and independent readback have now completed. Every field of
all 8,103,610 source hit rows and all 5,815,847 protein links matched the
original tables. All 526 taxon-level join summaries matched independent source
counts. There are 3,792,842 representative proteins with at least one GA hit
and 2,023,005 with no GA hit; non-detection remains explicitly qualified and
is not evidence of biological absence. Joining sequence-level hits back to
protein links produces 8,249,685 protein-linked hit rows. This differs from
the source hit count because identical sequences may have multiple protein
links and four marker-only queries lack representative links.

Completed receipts are archived in `metadata/full_domain_database_{receipt,readback}.json`.
The immutable database, source-shard manifest and checked 526-row
`taxon_domain_detection.tsv` remain under `results/domains/full-domain-database-v1/`.
The exported table hash and totals were checked again before archiving the
receipts. Overlap resolution and domain architecture analyses remain incomplete.


## Full overlap inventory (running)

`map_full_domain_overlaps.py` is now running across every query in the validated
database, under `fungal-full-domain-overlaps-20260917.service`. It emits all
pairs whose inclusive envelope intervals intersect, with their alignment
intersection sizes, same-family and same-clan flags, both directions of curated
NE relationships, both directions of coordinate containment, and their
conjunction. Pair flags remain separate: neither shared clan membership nor a
curated nested relationship automatically resolves an observed overlap.

The per-query inventory retains all source hits in alignment-coordinate order,
including model type and accession, explicit no-hit queries, counts of overlap
pairs and HMM coverage below0.70. Nonoverlapping raw annotations are labelled
as such, not as validated architectures. Overlap pairs are sorted and enumerated
with an interval sweep; no best-scoring hit is silently selected or discarded.
Output goes to `results/domains/full-domain-overlaps-v1/`.

`check_full_domain_overlap_sweep.py` compared 400 synthetic input sets to an
exhaustive integer-residue-set calculation under shuffled input order. All
32,997 overlapping pairs matched in identity, alignment/envelope size,
family/clan flags and directed curated-plus-coordinate nesting. This checks
the software kernel; independent readback of the full empirical output remains
required after completion. Architecture assignment, competing-model resolution
and evolutionary gain/loss interpretation are still pending.

The plan uses one CPU, a 4 GiB memory estimate and a 16 GiB service limit, with
50 GiB minimum free disk and a 20 GiB output allowance. The broad0.25–6 hour
planning window is not a measured ETA. Full configuration, pins, fixture
receipt and launch identity are versioned as `metadata/full_domain_overlap_*`.


The complete empirical overlap readback is now queued under
`fungal-full-domain-overlap-readback-20260917.service`.
`advance_full_domain_overlap_readback.py` waits on the exact producer PID and
creation time, requires its complete receipt and pinned plan, then runs
`readback_full_domain_overlaps.py`. The reader uses exhaustive within-query
combinations, independently of the interval-sweep kernel, to check the complete
pair universe and every alignment/envelope, clan, family, containment and NE
flag. It also checks all query inventories, ordered source annotations,
partial-HMM counts, status labels and aggregate totals. Extra, missing or
repeated pair/query rows fail the readback.

Three hundred additional synthetic cases (17,610 overlap pairs) agreed across
all nine pair annotations between the independent exhaustive calculation and
the production sweep. Full empirical checks remain pending. The reader uses
one CPU and a 16 GiB service cap (4 GiB planning allowance); its broad planning
window is0.25–8 hours after producer completion. No high-multiplicity query is
excluded. Config, fixture evidence and exact live service identity are in
`metadata/full_domain_overlap_readback_*`. A transient-service failure or
partial output requires review; source products are never overwritten.


The complete overlap inventory and exhaustive empirical readback have now
passed. Every query, ordered raw annotation and overlap pair was checked.
Both gzip artifact hashes were reverified before archiving the receipts in
`metadata/full_domain_overlap_{receipt,readback}.json`.

Across 5,713,603 sequence queries and 8,103,610 hits, there are 3,531,482
inclusive envelope-overlap pairs, of which 3,376,062 also overlap in alignment
coordinates. A total of 3,466,306 envelope pairs share a nonempty Pfam clan;
26,157 agree with both a directed curated NE relationship and coordinate
containment. These flags can overlap and are not mutually exclusive classes.
They do not establish which hits represent the biological architecture.

The query inventory retains 932,514 queries with overlapping hits requiring
review, 2,789,760 with nonoverlapping raw annotations and 1,991,329 without a
GA hit. These are sequence-query counts, not protein-link counts; they also
include four marker-only queries. Ordered annotations retain model types, so
Pfam Family, Repeat, Motif and other types are not silently relabelled Domains.
Choosing among competing models, validating nested segments, propagating
annotation uncertainty and comparing architectures on reconciled trees remain
required. The raw and checked overlap outputs are preserved unchanged.


## Clan competition implementation and planned architecture sensitivities

The [Pfam FAQ](https://pfam-docs.readthedocs.io/en/latest/pfam-faq.html)
describes selecting the lowest-E-value match among overlapping clan members.
The [2016 Pfam paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4702930/)
also distinguishes nested domains and a database-curation exception for small
overlaps that depends on family-wide frequencies. That exception is not a
blanket license to ignore every overlap shorter than 20 residues in this
project.

The implemented kernel `scripts/domain_clan_competition.py` uses within-query,
within-clan greedy competition, retaining a full disposition for every input
hit. The planned baseline uses alignment spans and independent domain
E-values; envelope spans and domain bit scores define three additional
sensitivity combinations. These are project policies, not a claim to reproduce
a particular PfamScan version. Decimal arithmetic preserves very small
serialized E-values. Secondary score and hit ID select a reproducible
representative, but ties on the primary ranking criterion remain explicit,
including rounded zero E-values. Cross-clan and missing-clan overlaps remain
unresolved. A chain of overlaps does not collapse to one component-wide winner.

A curated directed nesting relationship plus strict alignment containment
preserves both candidate hits. This is an annotation hypothesis: the kernel
does not infer physical nesting, remove internal sequence segments or recover
HMM match-state boundaries from a domain table. Equal boundaries do not meet
the containment exception. Model types, repeat instances and all suppressed
matches must remain available in the full output. No-hit queries remain
unknown for biological absence. Partial-HMM coverage and policy disagreement
must be propagated separately before phylogenetic architecture tests.

`check_domain_clan_competition.py` passed chain and input-order checks,
sub-floating-point E-value precision, rank ties, missing/different clans,
directed strict nesting, coordinate and ranking sensitivities, empty inputs
and invalid-input rejection. These are software checks. Full-data competition,
independent empirical validation, architecture calls and gain/loss analyses
have not yet run.
