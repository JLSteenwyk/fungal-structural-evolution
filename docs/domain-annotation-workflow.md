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
