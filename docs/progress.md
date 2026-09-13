# Progress and evidence

| Milestone | State | Completion evidence required |
|---|---|---|
| Local repository and original objective | Created | Git history and docs/objective.txt |
| GitHub remote | Published and public visibility verified via GitHub API | https://github.com/JLSteenwyk/fungal-structural-evolution |
| Catalog discovery | NCBI fungal catalogs inventoried; additional sources pending | metadata/source_receipts.json and metadata/catalog_summary.json |
| 500 fungi + 25 outgroups | 501 fungal entries +25 outgroups acquired; species identity review pending for 22 fungal labels | metadata/analysis_manifest.tsv; metadata/taxon_label_review.tsv |
| Data QC and structure inventory | Broad BUSCO and assembly-statistics collection complete for 526; exact-sequence structure inventory ongoing; lineage-specific and contamination QC pending | Per-taxon receipts; metadata/busco_dataset_receipt.json |
| Species tree and discordance | 125 profile and MAFFT alignments complete; 49,027-site profile and 63,750-site MAFFT matrices built; guide and gene trees running; support/discordance/sensitivity pending | metadata/initial_species_matrix_receipt.json; results/phylogeny/ |
| Families, domains and reconciliation | Representatives prepared for all 526; 64-taxon computational core complete; full assignment running; reconciliation validation pending; marker Pfam search and raw annotation complete; additional full-proteome search running | metadata/orthology_input_manifest.tsv; results/orthology/ |
| Structural atlas | Existing-model retrieval ongoing; 128 local ESMFold predictions independently verified; expanded queue running; direct comparisons executed for 865 marker taxon pairs and 738 domain comparisons; full atlas and clustering pending | metadata/direct_structural_comparison_receipt.json; data/structures/ |
| Evolutionary objectives 1–8 | Pending | Estimates, uncertainty, tests and sensitivity |
| Figures, methods and case studies | QC and exploratory geometry figures/methods available; evolutionary results and case studies pending | docs/figures/; docs/methods-draft.md |

## Chronological execution record

Earlier entries below describe the state at that time and are superseded by later evidence. In particular, early references to “unique species” counted distinct taxon IDs; the taxon-label review establishes that species-level identity is still pending for some entries.

Previous conversational turn supplied a goal prompt; it did not execute analyses. First execution turn initializes authoritative project records and obtains catalog evidence. Keep the full objective active.

Initial catalog evidence: 26,963 assembly records spanning 7,483 NCBI species IDs; 2,810 species have an annotation name in at least one record. This does not verify downloadable proteomes or balanced lineage coverage. Parser rejection checks passed; cached rerun preserved candidate and receipt checksums.

Taxonomy and sampling discovery progress: ranked-lineage snapshot joined to assemblies; coverage audit spans 19 named NCBI phyla plus unresolved entries. There are 35 missing species IDs requiring merged/deleted-ID checks. One provisional assembly selected for each of 2,810 annotated species, without committing a final sample. Imported 153 published timetree taxa with original sources and verified publisher archive checksum. Published sources must be checked for genome versus transcriptome provenance before outgroup selection.

Next: resolve missing taxonomy; validate candidate proteome availability and annotation metrics; select broad and transition-focused sets; identify approximately 25 actual genome-backed non-fungal outgroups. No structures or evolutionary estimates have yet been generated.

Current execution: corrected proteome availability checks running under tool session 42778 (must revalidate live handle on continuation). Initial process 54270 was explicitly terminated after detecting malformed URL generation. No failed malformed URLs count as biological missingness. Taxonomy audit resolves all 35 missing IDs as deleted; manual exclusion/replacement review remains. Added targeted close-outgroup assembly inventory, including independent Parvularia genome evidence.

Full-scale sampling/QC milestone: `metadata/fungal_sampling_draft.tsv` contains 500 unique fungal species from 18 NCBI phyla, 170 orders and 271 families. `config/sampling.json` and `scripts/select_fungal_candidates.py` reproduce selection. Ecological traits remain unknown, outgroups are still being curated, and Sanchytriomycota is missing pending external annotations.

Live acquisition handles at this update: 38081 (500-proteome downloads and checksum/FASTA validation), 42778 (2810-candidate HEAD availability checks). Revalidate these handles before restart. `scripts/summarize_downloads.py` records partial terminal results without marking unprocessed taxa failed. Protein counts are raw annotation records and include possible isoforms. No phylogenetic or structural analyses completed yet.

Availability scan completed: 2,806/2,810 current candidate URLs served successfully; four returned 404 (see proteome_availability.tsv). The original process's printed aggregate included historical malformed URLs; the output table correctly contains exactly 2,810 current URLs. GCA_900290405.1 (Saccharomyces jurei) is in the QC draft but has no protein file/checksum entry; find an alternative assembly or published annotation before final sampling.

Retrieved and validated both Sanchytriomycota proteomes (7,220 and 9,368 protein records) and their published genome files. Expanded draft includes 502 fungi, maintaining approximately-500 scope. Cross-source assembly/strain correspondence needs audit before CDS-based analysis. Initial Corallochytrium and Chromosphaera genome bundles retrieved with publisher checksums; additional external bundles remain downloading.

Handles to revalidate: 38081 (NCBI proteome download), 44022 (external genome bundles), 74981 (isolated BUSCO environment installation). Availability session 42778 is terminal and should not be restarted unless explicitly retrying updated sources. Prior goal turn was progress, with published data and active acquisitions.

BUSCO 6.1.0 installation verified via conda run; explicit Linux package lock committed. Downloaded and hashed eukaryota_odb12.2 (125 markers, dataset version 01, creation 2026-05-13). Prepared 408 currently available fungal proteomes for initial raw-proteome QC; excluded internal-stop/noncanonical records are counted in qc_input_receipts.json, terminal markers removed with audit. This is a rolling full-scale batch, not a separate pilot. Remaining taxa will enter the same QC workflow when downloads finish. Lineage-specific scores and isoform-aware assessment remain pending.

External bundles completed after allowing audit of terminal dot markers in Abeoforma (6,131 terminal-dot records; no internal dots). Original files remain unchanged and require normalization for analysis. Creolimax and the four other unicellular holozoan genome bundles now have publisher-checksum evidence. No final 25-outgroup manifest yet.

Saccharomyces jurei annotation recovery: cited ENA ERZ491603 endpoint returns 404; taxon-specific UniProt proteome search returned no entries; NCBI assembly GenBank file contains zero CDS and zero translations. Original missing proteome is not a transient HEAD-only issue. Keep this taxon unresolved while searching published annotations; do not silently replace it with another species.

Latest verified execution state: original 500-proteome acquisition is terminal (499 validated, 1 unresolved S. jurei). External acquisition and BUSCO installation are terminal success. BUSCO batch handle 94351 is live; first two taxa returned zero exit codes and produced summary outputs. Do not restart the same batch while that handle is live. Once terminal, prepare all available 501 fungal inputs (499 NCBI plus 2 gap-fill taxa) and rerun to process missing taxa, reusing successful outputs. Original 502 draft includes unresolved S. jurei.

Outgroup milestone: config/outgroups.json defines 25 genome-backed candidate species: one nucleariid, nine unicellular holozoans (four NCBI plus five published bundles), one apusomonad, five amoebozoans and nine animals. All 20 NCBI proteomes passed checksum/FASTA validation; the five external bundles are verified, with terminal-marker normalization explicitly required for Abeoforma. Closest-relative representation can improve if additional nucleariid annotation sources become usable; do not substitute fungal boundary taxa for non-fungal outgroup count.

Combined sampling_manifest.tsv now records 502 fungal candidates and 25 outgroups (527 total), with explicit provisional/QC states and unresolved S. jurei. It is not frozen. Shared eukaryote BUSCO summary table contains completed jobs only with marker-count consistency checks. Initial batch 94351 remains live; preparation process 70767 adds remaining taxa while reusing unchanged input files. Never launch overlapping BUSCO work on the same taxon.

Working analysis set is now explicitly 501 available fungi + 25 outgroups (526 total), within the approximate target. S. jurei is retained in the 527-candidate manifest but excluded from analysis_manifest.tsv with evidence in sampling_exclusions.json; recovery can add it later. This does not block the full-scale analysis.

AFDB inventory completed (session 39260 terminal): 429 candidate taxa have legacy bulk archives at queried species/assembly IDs, 97 have none at those IDs, and one lacks a resolved taxonomy ID. Both v3 and v4 suffix versions are present; selected_archives deduplicates shard versions. These counts do not establish exact reuse or absence from current AFDB. BUSCO initial batch 94351 remains live; 104 successful summaries were inspected, with completion 74.4–100% in that partial, taxonomically ordered subset (not a whole-dataset quality conclusion). Remaining 118 prepared taxa await the next non-overlapping batch.

AFDB bulk retrieval attempt: session 34508 was explicitly terminated after verifying 401 responses from listed media URLs and 403 AccessDenied from the documented public object URL, with and without generation pinning. No structural archive was downloaded. Bucket metadata listing is accessible but anonymous object retrieval is not; this is an access limitation, not absence of structures. Preserve archive plan and use current services instead.

Fallback exact sequence reconciliation started under handle 88864: UniProt taxonomy-scoped sequence tables are matched exactly to prepared input sequences; AlphaFoldDB crossrefs nominate candidate structures. Current AFDB prediction API responded successfully in an independent read-only check. Structure sequence/coverage/confidence still need verification before reuse. BUSCO handle 94351 remains live.

Actual structural atlas acquisition began: current AFDB API and CIF files are accessible. Retrieval handle 94665 processes the first approximately 49,900 accessions available from completed UniProt matches, while matching handle 88864 continues across all taxa. Each accepted monomer requires exact API/CIF/input sequence equality and complete CA residue coverage. Per-model pLDDT summaries, version, source, checksums and protein mappings are retained. Unmatched/fragment models are recorded separately; their regions may be reused later after explicit coverage mapping. PAE is not yet downloaded, so domain orientation conclusions are not justified.

BUSCO handle 94351 is still live; additional prepared taxa await the non-overlapping continuation. Existing coordinate retrieval does not replace orthology, ecological curation, new structure prediction or the eight evolutionary analyses.

Phylogenetic marker extraction and bounded MAFFT alignment are now implemented. A staging-only extraction from 296 completed taxa verified 34,722 full-length source-identical sequences across all 125 eukaryotic markers; output is `results/phylogeny/markers-staging-20260913-a`. Downstream alignment rejects this incomplete snapshot. Three integrity tests pass, covering missing-taxon gating, changed source sequences, and alignment identity/residue preservation. No phylogeny has been inferred. The documented full-panel workflow awaits completion of all 526 BUSCO jobs, followed by alignment assessment and model/sensitivity decisions.

Latest committed snapshots report 336 successful BUSCO jobs and 807 verified AFDB models mapping to 814 input proteins in one taxon. This first-taxon concentration reflects retrieval order and is not a cross-fungal coverage estimate. Live handles remain 94351 (initial 408-taxon QC batch), 88864 (UniProt matching), and 94665 (current AFDB coordinates). The remaining 118 prepared QC inputs still require a subsequent non-overlapping batch.

Initial BUSCO session 94351 is terminal with exit code zero and all 408 successful receipts. Continuation session 43450 now runs the remaining 118 prepared inputs, reusing the existing 408 results. A batch lock now prevents future overlapping invocations. The full 526-taxon marker extraction remains pending its completion.

Assembly-matched NCBI GFF acquisition is running as session 10525, two workers across 519 taxa; the seven published-bundle taxa need separate annotation reconciliation. At the latest inspection 175 GFFs had passed publisher checksum and feature validation. Protein-to-gene mapping completed on a snapshot of 61 taxa: all 632,728 input proteins mapped to exactly one gene, representing 629,835 genes; 2,657 genes have multiple protein products. This is not a duplication or alternative-splicing result. Isoform selection has not yet occurred. Mapping tables retain accession, local gene IDs, partial-CDS flags and source checksums. Six integrity tests pass, including multipart CDS, multi-parent ambiguity and cycle rejection.

Current committed coordinate snapshot: 1,300 exact-sequence verified models linked to 1,307 input proteins in one taxon. UniProt matching session 88864 and AFDB session 94665 remain live. Gene mapping and acquisition scripts, the annotated workflow, and resource allowance are now documented; phylogenetic and evolutionary inference remain outstanding.

Representative-protein selection is implemented with a deterministic longest-per-uniquely-mapped-gene baseline and lexical accession ties. The first staging run (61 taxa) retained 629,835 proteins and recorded 2,893 alternative products without deleting source sequences. Results and selection receipts are under `results/gene_representatives/staging-v1`; no orthogroup inference uses this incomplete staging set. Unresolved gene mappings are retained independently and excluded from direct gene-copy interpretation. Eight integrity tests pass, including distinct-locus preservation and published transcript-ID formats.

Mapping now supports the four published outgroup GFFs via exact mRNA IDs and Creolimax via CDS transcript/gene attributes. Updated mapping snapshot session 30163 is still running; do not launch another writer before it finishes. Two sanchytrid protein headers contain genomic coordinates matching the supplied genome contig naming, but correspondence and translation require verification before these can support locus assignments. Latest live inspections: BUSCO continuation 43450 at 468 successful/reused jobs; annotation acquisition 10525 at 342 validated files; UniProt 88864 and coordinate acquisition 94665 continue. Previous turn made concrete progress through executed representative selection and annotation-parser integration; overall goal remains incomplete.

NCBI GFF acquisition session 10525 completed successfully: all 519 annotations validated, totaling 1,183,724,940 compressed bytes. Receipts are versioned in metadata/annotation_download_receipts.json. Mapping session 30163 is terminal after processing a 304-taxon snapshot. Four external outgroup GFFs map all input proteins; Creolimax retains 50 unresolved protein IDs. Inspection shows these proteins use gene IDs whereas the primary GTF mapping uses transcript IDs; this source-specific naming case needs resolution rather than accession truncation.

The sanchytrid coordinate audit executed across all 16,588 proteins: 7,220/7,220 Amoeboradix and 9,367/9,368 Sanchytrium translate exactly from the deposited genomes. Node/length prefix aliases resolved coverage-suffix differences only when translation matched. One Sanchytrium protein (NODE_1518_length_898_cov_0.933667.p2) specifies positions 896–510 on an actual 837-base contig and remains unresolved. Verified CDS FASTAs and mappings have source/output checksums; provisional ORFs are not equated with resolved gene models. Ten integrity tests pass, including coordinate bounds and reverse-strand translation.

Full mapping/representative workflow session 66014 started after session 30163 exited; it includes all available annotations and the provisional sanchytrid ORF mappings. Its outputs are logs/gene_mapping_full.log, logs/gene_representatives_full.log and results/gene_representatives/full-v1. Do not start a competing writer. Marker workflow session 24260 is queued behind the live BUSCO batch lock, then will require complete successful QC before extraction and bounded MAFFT alignment. At latest inspection 525/526 BUSCO jobs had succeeded; Branchiostoma floridae (O2700040) remained running in batch 43450. UniProt matching and AFDB retrieval remain active under their existing handles. Overall evolutionary analyses remain pending.

Full BUSCO batch 43450 completed with exit zero: all 526 taxa now have successful QC receipts. Marker workflow 24260 extracted 59,840 source-identical single-copy sequences across 125 markers and launched bounded MAFFT alignment. The first six alignments completed and validated; the remaining alignments are still running. Full mapping/representative workflow 66014 also remains live.

The full-dataset QC summary and figure are generated in results/qc and displayed from docs/figures/marker_recovery.svg. All 526 taxa are represented. Microsporidia show median 42.4% complete-marker recovery; two published outgroups (Pirum, Abeoforma) each recover 5.6%. These require marker/lineage and taxon sensitivity analyses; no blanket exclusion was made. The figure was visually inspected. The methods draft now distinguishes executed preparation from pending inference.

Mapping-gap audit: Phaffia rhodozyma F264483 has 6,380 proteins and no gene features in its NCBI GFF (CDS/mRNA records are present); its unresolved mappings are not evidence of missing proteins. Creolimax's 50 extra products use exact gene IDs. The parser now supports those gene IDs, with a --taxon option for a focused refresh. After live workflow 66014 is terminal, refresh OFS1403592, then create an updated representative snapshot (do not overwrite full-v1). This correction has passing format tests but has not yet been applied to the running older parser process. Ten integrity tests pass overall.

Alignment assessment is implemented and executed. Initial 11-marker MAFFT audit retained a median 7.4% of columns at 50% occupancy, reflecting long insertion-rich full-protein alignments. This prompted a full-panel profile-based comparison, not a separate pilot. HMMER 3.4 profile alignment session 76822 completed successfully across all 125 markers, validating exact source sequence preservation in full Stockholm files and profile-length consistency in match-state outputs. Full profile audit session 15516 is terminal: 50,941 match-state columns, 49,027 retained under the candidate 50% occupancy mask (median per-marker retained fraction 98.27%). No filtering policy has yet been used to infer a tree.

Profile outputs are results/phylogeny/profile-alignments-full-v1 and profile-alignment-audit-full-v1. Audit statistics and receipt are copied to metadata/profile_alignment_statistics.tsv and metadata/profile_alignment_audit_receipt.json; refresh those copies after rerunning the documented audit command. All 11 integrity tests pass, including informative-column and ambiguity handling. MAFFT workflow 24260 remains live for alignment sensitivity; full gene mapping/representative workflow 66014 remains live. Structure retrieval and UniProt matching continue under their existing handles. Next phylogenetic step: construct audited concatenated matrices, assess taxon coverage and fit trees with model and alignment sensitivity.

Full mapping phase of 66014 has now finished for all 526 taxa; representative selection in the same session is still running. Snapshot totals: 5,927,745 proteins, 5,904,721 unique-gene mappings, 16,587 provisional ORFs, 6,436 unmapped and one multi-gene mapping. There are 5,792,823 uniquely mapped genes, of which 59,504 have multiple protein products. The one multi-gene protein is in Capitella teleta O400682. Phaffia accounts for 6,380 unmapped proteins; Creolimax's pending naming correction accounts for 50. These are mapping classes, not evolutionary gene-duplication estimates.

Full representative session 66014 completed successfully. Creolimax naming was then refreshed with `map_proteins_to_genes.py --taxon OFS1403592`; corrected representative session 8824 also completed, writing full-v2. The final current baseline retains 5,815,847 proteins across 526 taxa, with 111,898 alternative products retained in source records. Mapping and representative receipts are versioned. This resolves the 50 Creolimax naming gaps without changing representative count; other unresolved/provisional classes remain explicit.

The full profile matrix was built and validated: 526 taxa, 125 markers, 49,027 columns, with missing-data padding, taxon coverage, site coordinates and partition records. An initial unpartitioned LG+F+G4 guide tree is running as session 98857, 16 threads and 32 GB limit; IQ-TREE reports approximately 15.9 GB required. This is not a final supported tree. Composition and low-coverage concerns are documented in the phylogenetic workflow. All 12 integrity tests pass, including concatenation coordinate/missingness checks.

OrthoFinder 3.1.5 official archive was downloaded and checksum-verified, then installed in an isolated Python 3.12 environment (installation session 23453 terminal). Help invocation succeeds; dependency locks and software receipt are versioned. Core/assignment workflow documentation was reviewed from official sources. Orthology core selection and actual inference have not yet been launched. MAFFT session 24260 and initial guide 98857 remain running; representative preparation is terminal. Overall evolutionary objectives remain incomplete.

Orthology inputs were frozen for all 526 taxa. The 64-taxon reference core contains 54 fungi and 10 outgroups, 597,213 proteins, all 19 represented fungal groups and eight outgroup categories. The 462 additional taxa contain 5,218,634 proteins. Source checksums and stage labels are in orthology_input_manifest.tsv; no taxa were dropped. Selection was checked for input-order independence and agreement with the frozen manifest. This reference/assignment split is a computational dependency of the full study.

Core inference is running as session 62364 (`run_orthology_core.py`), with 32 search threads and eight analysis workers; required runtime checks for MCL, FAMSA and FastTree passed, and search preparation began. Results are under results/orthology/core-v1/Results_Sep12; control/config/logs under results/orthology/core-control-v1. Resources were rechecked before launch (about 948 GB RAM available, 12 TB disk free). Full assignment requires successful core output and an updated resource estimate; it has not yet begun.

Other jobs revalidated live this turn: initial guide 98857, MAFFT 24260, UniProt matching 88864, and AFDB coordinate acquisition 94665. Latest inspected logs show 124 taxa processed by matching and over 5,500 verified model retrievals; these are running progress counts, not completed atlas coverage. No structural evolutionary results are yet inferred.

Individual marker-tree inference is implemented for all 125 profile markers, including explicit per-marker coverage filtering, model comparison and SH-aLRT support. First launch 37936 is terminal: IQ-TREE rejected an option spelling and returned errors before inference. Corrected run 61627 writes marker-gene-trees-v2, uses documented single-dash frequency/rate options and has progressed through likelihood optimization. Four jobs run concurrently at two threads each; resource allowance is documented. All 13 integrity tests pass, including the new coverage-rule test. No individual gene tree is yet claimed complete.

Reviewed official AlphaFold bulk instructions: a Google Cloud account is required for public proteome archive downloads, while those dataset downloads are described as no-cost. The current PATH has neither gcloud nor gsutil. An optional asynchronous question asks whether the user can authenticate an existing account; this does not pause current individual-file retrieval or other project work. No paid resources or accounts were provisioned. This is a possible acceleration path, not a whole-project blocker.

Structure-to-phylogeny mapping executed across the full 59,840-marker-protein set. Snapshot results/structural_markers/snapshot-v1 links 151 proteins from two taxa to 57,413 species-matrix positions, retaining coordinate and pLDDT provenance. None of the markers yet had four structurally represented taxa; no branch result is claimed. Matrix/profile/CIF sequence and checksum checks validate these mappings. Fifteen integrity tests pass, including insertion/gap residue indexing and taxon-balanced queue behavior.

To broaden comparative coverage sooner, retrieval now prioritizes marker accessions round-robin across taxa before all remaining atlas candidates. Old session 94665/PID 936663 was intentionally terminated for this scheduling change, confirmed terminal with exit 143, and all 6,924 retained JSON receipt records parsed successfully. Replacement session 73872 revalidated cached model checksums and started a 1,043,924-accession snapshot, including 12,007 priority marker accessions across 113 taxa; 183 taxa had completed matching at queue creation. It retains the same two-worker request concurrency and adds bounded futures plus clean-stop handling. Matching 88864 continues; no paid resources or authentication changes occurred.

The updated structural mapping (session 14466 terminal, snapshot-v2) verifies the scheduling benefit: 452 marker proteins, 118 taxa, 422 distinct models and 159,837 matrix-residue links. Fifty-four markers now have at least four linked taxa, compared with none in snapshot-v1. The current metadata copies refer to snapshot-v2; both raw snapshots are retained. Confidence/domain/source assessment and supported phylogenies are still needed before branch inference. Measured short-window retrieval throughput and coordinate storage estimates were added to the resource record.

Direct structure/sequence comparison session 57235 completed successfully: 2,209 qualified rows across pLDDT 50/70/90, representing 865 within-marker taxon pairs across 111 markers; 386 threshold-specific rows were excluded by coverage. At threshold 70, 858 comparisons passed. Proper-rotation RMSD and explicitly defined local Cα distance-change metrics use matched phylogenetic positions; no branch-rate or independence assumption is made. Matrix provenance was independently rechecked during execution, and the reusable script now enforces it directly. All 16 integrity tests pass, including rigid-rotation/translation invariance and reflection exclusion.

The direct-comparison figure was generated and visually inspected; artifacts are results/structural_comparisons/snapshot-v1, with an SVG copy under docs/figures. The plot is descriptive only and includes no correlation significance test. Global-versus-local discrepancies motivate PAE, domain-orientation and alignment review. Supported gene/species trees, core orthology and structural retrieval remain in progress; no structural evolutionary rate result is claimed.

PAE acquisition and confidence sensitivity completed: all 422 distinct models in structural mapping snapshot-v2 have verified, version-matched PAE matrices (54,375,930 compressed bytes; zero failures). All 858 qualifying pLDDT ≥70 marker taxon pairs were assessed at PAE thresholds 5, 10 and 15 Å, requiring both directions in both models, for both all-nonadjacent and local-neighborhood residue pairs. This produced 5,148 sensitivity rows with zero PAE exclusions. The analysis independently reproduced every baseline RMSD and local residue-pair count before filtering. Receipts: `metadata/marker_pae_retrieval_receipt.json`, `metadata/pae_sensitivity_receipt.json` and `metadata/pae_sensitivity_figure_receipt.json`. Full results and all-pair manual review ordering remain under `results/structural_pae/comparisons-v1`; SVG copied to `docs/figures/pae_sensitivity.svg` after visual inspection. All 18 integrity tests passed, including directional PAE and residue-index checks.

At PAE ≤10 Å, the median retained fraction of nonadjacent residue pairs is 0.7455. The across-comparison median of mean absolute distance changes is 0.4849 Å before filtering and 0.3064 Å among confident pairs. These averages use different residue-pair compositions and are not tests of an evolutionary effect. For the highest-global-RMSD comparison (marker 5000823at2759, the explicitly labeled Saccharomyces hybrid versus Pachysolen tannophilus), global RMSD is 34.7439 Å; only 33.32% of nonadjacent residue pairs pass PAE ≤10 Å, and mean absolute distance change is 8.4510 Å unfiltered versus 1.0132 Å on confident pairs. This motivates domain/orientation and annotation review; it does not establish an artifact, structural acceleration or a hybridization effect.

Taxon identity correction: all 526 labels were screened and 22 fungal entries were flagged (21 incompletely identified labels and one explicit hybrid). The earlier description of 501 unique fungal species was too strong: the working set contains 501 fungal entries, and independent species-level identity remains under review. `metadata/taxon_label_review.tsv` records the flags and required sensitivities. Existing data and running inferences were preserved. Initial guide-tree inference, marker-tree inference, MAFFT alignments, full-input OrthoFinder core inference, UniProt matching and AFDB coordinate acquisition were each confirmed live during this work; none is presented as completed evolutionary inference.

The complete 125-marker MAFFT batch (session 24260) is now terminal with exit code zero. An audit attempted before the final marker finished correctly rejected the incomplete batch; after verified completion, the full audit completed successfully. It records 697,866 raw full-protein columns and 63,750 columns at 50% occupancy. The alternative MAFFT concatenation completed for all 526 taxa, preserving explicit missingness, 125 partitions and residue-to-matrix coordinate provenance. Receipts and statistics are copied to `metadata/mafft_alignment_audit_receipt.json`, `metadata/mafft_alignment_statistics.tsv` and `metadata/mafft_species_matrix_receipt.json`. The profile matrix remains 49,027 sites; no new tree inference has been launched from the alternative matrix yet.

Residue-level alignment correspondence completed for all 125 markers. The 50%-occupancy masks retain 22,778,444 profile residues and 28,568,763 MAFFT residues, with 22,269,691 identical protein residues in common. Among this common residue universe, the profile alignments define 5,113,546,484 cross-taxon residue-pair edges, MAFFT defines 5,021,944,588, and 4,871,756,694 edges agree. Pooled edge Jaccard is 0.925532; the median across markers is 0.944903 (range 0.573155–0.993988). Marker 530740at2759 has the lowest correspondence and is a priority for alignment sensitivity, without assuming either method is correct. The method conditions on shared retained coverage, so these values are not whole-protein alignment accuracy. Complete marker statistics are in `metadata/alignment_correspondence.tsv`, hashes in its receipt, and per-column diagnostics in `results/phylogeny/alignment-correspondence-v1/profile_column_correspondence.tsv`. The SVG was visually inspected and copied to `docs/figures/alignment_correspondence.svg`. All 20 integrity tests passed, including exact checks of the edge-counting algorithm against explicit enumeration. Species-tree topology/support sensitivity and alignment effects on structural comparisons remain outstanding.

Pfam domain annotation preparation completed and production searches started. Release 38.2 (30,134 profiles; UniProtKB 2025_03 basis) was acquired from the version-specific EBI archive with pinned publisher MD5s and local SHA256 checksums. Five compressed resources total 419,509,607 bytes and decompress to 2,253,592,241 bytes before the duplicate HMM chunk storage. All profiles have unique accessions, lengths and gathering thresholds; all 30,134 metadata entries validated. Types comprise 15,887 Domain, 12,074 Family, 1,560 Repeat, 314 Coiled-coil, 166 Disordered and 133 Motif entries. `metadata/pfam_release_receipt.json` records sources, checksums and 64 deterministic search chunks.

All 59,840 full marker-protein records across 125 markers and 526 taxa were reconciled against their source sequences and collapsed to 58,883 exact unique sequences (33,047,630 residues) for search only. All taxon/protein associations are retained. `metadata/marker_domain_input_receipt.json` pins this input snapshot. Session 25652 runs `run_marker_domains.py` at `results/domains/marker-search-v1`, with four HMMER 3.4 `hmmsearch --cut_ga` processes, each using two HMMER worker threads. All four child processes were observed actively using CPU. The binary, configuration and sources are pinned in `metadata/marker_domain_search_config.json`; raw outputs and per-chunk receipts remain outside Git. Neither the domain search nor interpreted annotation is complete yet.

A complete-search summarizer is implemented and tested to preserve Pfam types, versioned identities, alignment/envelope coordinates, HMM coverage, gathering scores and overlapping hits. It records undetected profiles without claiming domain absence and does not resolve architectures automatically. All 22 integrity tests passed, including HMMER table orientation, inclusive coordinates and overlap boundaries. A format check validated 100 complete rows from the active production output against real Pfam metadata and target lengths; it is not an annotation snapshot or a biological pilot. Full-proteome domain annotation, overlap/architecture resolution, domain-restricted structural comparisons and evolutionary domain analyses remain outstanding.

Assembly-quality collection completed for all 526 taxa: all 519 assembly-version-matched NCBI statistics reports passed publisher checksum, accession and format checks (3,834,008 bytes; no failures), and all seven external genome FASTAs passed checksum/DNA validation and direct record-level measurement. The NCBI table preserves whole-assembly and primary-assembly scopes separately. Whole-assembly contig N50 is reported for 474 assemblies; the 45 missing fungal values remain unknown, not zero. The seven external FASTA-record N50s are not assumed equivalent to NCBI contig N50. NCBI report retrieval completed under session 51258; a subsequent hash-verified cache-only refresh clarified that external sources are documented in their separate metrics table. External measurement session 92095 completed successfully.

The combined quality table includes all 526 taxa and is linked to source-verified BUSCO summaries. Its figure has 474 NCBI points and seven external points in separate panels and was visually inspected. Tables and receipts: `metadata/assembly_quality_metrics.tsv`, `metadata/external_assembly_quality.tsv`, their receipts and `metadata/assembly_quality_figure_receipt.json`. Figure: `docs/figures/assembly_quality.svg`; full outputs: `results/qc/assembly-quality-v1`. Pirum and Abeoforma combine record N50s of 1,637/2,170 bp with 5.6% broad marker completeness, strengthening the need for detection-aware and taxon-exclusion sensitivities before interpreting absences. No taxa were automatically excluded and no contamination-free status was assigned. All 24 integrity tests passed.

Ecological evidence preparation completed for the current snapshot. The pinned publisher FungalTraits workbook gives exact, unambiguous genus-level candidates for 482 fungal entries; 19 have no exact genus match and the 25 outgroups are unassigned. All 526 candidate rows retain species verification and confirmatory eligibility as false. Primary/secondary lifestyles and trait strings remain separate, missing traits remain unknown, and main-manifest ecology fields were not overwritten. The source workbook contains 10,771 genus rows and six duplicate names; duplicate matches are not resolved automatically. `metadata/ecology_candidates.tsv` and its receipt preserve source rank and coordinates.

Six primary-literature statements were then curated separately in `config/ecology_species_evidence.json` and joined to selected taxa in `metadata/species_ecology_evidence.tsv`. The species-level evidence flags two conflicts with genus-based ECM candidates: Amanita inopinata and A. thiersii are classified as asymbiotic by Hess et al. 2018. The A. muscaria source is variety-specific and remains flagged for taxonomic-concept linkage; selected-isolate experimental validation is not claimed. Amanita transition labels preserve the source's single-origin interpretation instead of counting related species as replicate transitions. These statements do not complete ecological coverage or association testing. All 25 integrity tests passed; source DOI/page locators and limitations are documented in the ecological evidence workflow.

The first three completed marker trees (115730at2759, 129234at2759 and 201273at2759) passed a new support/provenance audit. They contain 481, 472 and 478 taxa and 235, 166 and 254 aligned columns, respectively. Across 1,422 internal branches, 1,390 have reported SH-aLRT values and 32 are unreported; missing values remain distinct from zero. Reported-support medians are 84.9, 82.2 and 84.55. Raw trees retain all taxa and branches. Canonical split records include each tree's full taxon universe and are invariant to arbitrary unrooted display orientation. This is an explicitly incomplete 3/125 snapshot, with 122 pending markers listed in its receipt; it is not a species-tree or discordance result. All 27 integrity tests passed. Summary: `metadata/marker_tree_support_snapshot.tsv`; detailed data: `results/phylogeny/marker-support-snapshot-v1`.

Full representative-proteome domain input preparation completed under session 93523. All 526 source FASTAs passed their representative-receipt checksums and count checks; all 5,815,847 protein/taxon records are retained in the link table. Exact sequence deduplication yields 5,713,599 unique full-proteome sequences: 58,879 already occur in the marker input and 5,654,720 additional sequences contain 2,520,735,813 residues. The marker partition serves 60,033 full-proteome protein records; its hits are reusable only after the complete marker search is validated. Four sequences in the marker input are absent from the representative baseline, so marker and representative input universes are not assumed identical. No gene-copy or species identities were collapsed by this search preparation. `metadata/full_domain_input_receipt.json` pins the inputs and output FASTA/link tables under `data/domains/full-inputs-v1`.

A resource projection from 45 completed Pfam production chunks gives 44.57 hours at current concurrency under linear residue-count scaling, with an explicitly provisional 2–7 day planning range. `metadata/full_domain_search_estimate.json` records the source evidence. Additional full-proteome searches are not yet launched; domain architecture interpretation and evolutionary analyses remain pending. Raw E-values will retain their search-partition context when marker and additional results are joined.

Independent read-back verification of the completed full-domain inputs passed under session 2694: all additional FASTA identifiers matched their actual sequence SHA256 values; identifiers were unique and disjoint from the marker-query set; every one of the 5,815,847 protein links resolved to its declared query partition; and all per-taxon counts, residue totals and marker-reuse counts matched the completion receipt. This validates the prepared input catalogue, not downstream annotations.

Additional full-proteome domain execution is now queued under session 39822, waiting on the live marker-search lock. The dependency chain first runs the complete marker annotation summarizer and launches `run_full_domains.py` only after that succeeds. The new runner verifies all marker profile chunks, result hashes, annotation artifacts, full-input linkage and identical Pfam/HMMER versions before starting four two-worker-thread jobs against 5,654,720 additional unique sequences. It reuses the marker CPU allocation rather than adding another concurrent domain batch. `metadata/full_domain_execution_plan.json` records the submission, source hashes and the existing 2–7 day resource estimate. Queue submission is not evidence that the full search has begun or completed. Its log is `logs/full_domain_search_queue.log`; eventual search outputs are `results/domains/full-search-v1`.

Prediction environment inspection: the existing `esmfold2` environment has Python 3.12, PyTorch 2.6.0+cu126 and Transformers 4.57.6; CUDA availability was verified without loading a folding model. Its `esm` package is the newer package, not evidence that FAIR ESMFold v1 is installed. The separate `openfold` environment resolves a local OpenFold source checkout but has no `esm` or Transformers module. The inspected default Torch cache contains ESM-2 650M weights, which do not establish availability of a complete folding model. A neighboring project has an ESMFold2-Fast loading example, but this project's model weights, inference compatibility, provenance and resource envelope remain unverified. No fungal structures were predicted during this inspection. GPU 1 was idle with approximately 49 GB total memory; GPU 0 had an active unrelated workload. Availability must be rechecked before inference. Marker Pfam search had 62 of 64 chunk receipts at this checkpoint; the validated full-proteome search remains queued behind completion.

Marker Pfam search finished all 64 chunks and all 30,134 models; complete annotation validation passed. There are 163,650 raw hits across 58,354 of 58,883 unique marker sequences and 13,983 sequences with overlapping hits. Types: 89,273 Domain, 30,904 Family, 42,111 Repeat, 921 Coiled-coil, 435 Disordered, 6 Motif. Raw hits are not resolved domain architectures. The queued full-domain runner passed its compatibility gates and is executing four HMMER jobs with two worker threads each against 5,654,720 additional unique sequences. Completion of marker annotation does not establish full-proteome search completion.

Local structure prediction preparation now verifies cached ESMFold v1 safetensors against official revision ba837a39b67e59941c3f017d6c2a064f567038d9 (SHA256 9a865162cdcaac8d5385c908fd0c620fccf0405adb3a5b49566a8d88831977ac). Small model files passed publisher Git-blob checks. A frozen marker queue from 363 completed taxon inventories retains all 59,840 taxon/marker links: 20,480 unique sequences without a reuse candidate, 13,170 nominated retrievals pending, 8,464 with verified reuse receipts, 16,769 awaiting inventory. The first group includes 10,536 sequences no longer than 512 residues, of which 14 contain noncanonical residues; 10,522 are eligible for the first local configuration and 9,958 are explicitly deferred by length/alphabet. A 128-sequence production chunk was launched on the independently checked idle GPU 1; six artifact-validator tests passed, including detection of confidence scaling errors, sequence/numbering errors and invalid PAE. Model loading and actual inference still require runtime verification. The short-protein chunk is part of the full run, not a separate biological pilot; it does not complete the atlas. See docs/structure-prediction-workflow.md for the pre-launch resource envelope and settings.

Runtime verification supersedes the launch-only prediction entry: all 128 first-chunk predictions completed and passed independent Biopython/PDB and NumPy artifact readback, preserving 129 taxon-marker links in 129 taxa and 15 markers. Lengths 70–255 residues, median inference 0.944 seconds, maximum allocated GPU memory 8.92 GB. A separate checkpoint-loading audit confirmed only two missing contact-regression parameters; a rejecting forward hook proved the contact head was not invoked during a complete repeated prediction. Repeated PAE and CA confidence matched exactly for that sequence. This is computational repeatability, not experimental accuracy. The identical configuration resumed for all 10,394 remaining canonical marker candidates no longer than 512 residues (session 61234); 9,958 other candidates remain length/alphabet deferred. The first chunk receipt is preserved at results/predictions/chunk-receipts/esmfold-marker-v1-chunk1.json.

OrthoFinder core finished successfully with all 64 expected tree taxa, 597,213 genes, 510,318 genes assigned to 24,215 orthogroups, and 86,895 unassigned. The complete 462-additional-taxon assignment controller was launched (session 1748) after input/software/core validation, with 32 search threads and eight analysis workers, DIAMOND, FAMSA and FastTree. Outputs are named full526v1 beside the core. Exact final 526-taxon scope and protected core artifacts are checked on computational completion; full scientific orthology/reconciliation validation remains pending. Pre-launch resource assumptions are documented in docs/orthology-workflow.md.

Assignment launch correction: session 1748 is terminal with status failed_or_requires_output_review. OrthoFinder 3.1.5 called removed `numpy.chararray` under NumPy 2.5.3 during reference-profile preparation and printed a traceback despite returning zero. No full assignment result was accepted. NumPy 2.2.6 is now pinned only in the isolated OrthoFinder environment, pip check passed, and the previously failing matrix function round-tripped a real alignment. The original core dependency lock and failed run are preserved. The revised controller validates the NumPy version and records the dependency-lock hash, and now rejects incomplete results even if OrthoFinder returns zero. A new full526v2 assignment run has started; output and configuration are in assignment-control-v2.

Domain/structure integration completed for the immutable AFDB snapshot: all 858 pLDDT ≥70 direct-comparison baselines were reproduced, yielding 738 conserved-domain comparisons across 463 marker–taxon pairs and 84 Pfam domains. Annotation eligibility retains single-instance, non-overlapping Pfam Domain hits with >=50% HMM coverage; coordinate comparisons require >=30 qualified matched sites and >=50% shared domain coverage. There are 381 marker-pair no-shared-eligible-domain exclusions and 95 domain coverage exclusions. Three mechanical geometry tests passed, and the plotted figure was visually inspected. Median identical-domain-site RMSD is 0.935 Å under the whole-marker fit and 0.509 Å under an independent domain fit; fitting improvement is mathematically expected and is not itself proof of movement.

Interdomain PAE/distance comparisons completed for all 375 eligible domain-pair/taxon-pair combinations (188 marker–taxon pairs), evaluated at 5/10/15 Å in both directions in both models. At 10 Å the median retained interdomain fraction is 0.3513 and 33 combinations have no passing pairs. A PH_SPT16 example (Pseudovirgaria hyperparasitica vs Pachysolen tannophilus, marker 4975094at2759) has 34.718 Å RMSD on domain sites under the global fit versus 0.577 Å under the domain fit. Its local within-domain pairs all pass PAE10, but none of 19,240 pairs to the Peptidase_M24-like region pass. It is recorded as an artifact-control candidate, not an evolutionary discovery. See docs/domain-structure-comparisons.md and tracked source/figure/placement receipts. Full species/gene phylogenetic inference, expanded structural atlas and all eight evolutionary analyses remain incomplete.

Native structural-alphabet extraction and feature auditing completed for all 422 models of the AFDB comparison snapshot (215,518 residues). Foldseek is pinned to commit e3fadcd07f971e864c094ac4f3a78bf4ed845e07 with executable/source hashes. Native amino-acid exports, 3Di state strings and ten-feature descriptors agree for every model; independent scripted re-extraction reproduced all 422 records exactly by identity. The audit reconstructed virtual centers, partners and geometric features from original full-backbone coordinates and matched native four-significant-digit descriptors. All 844 terminal states are explicitly invalid even though the native encoder assigns them an ordinary coil-state symbol. Counts passing focal pLDDT70, six-residue pLDDT70, and six-residue pLDDT70 plus maximum-context PAE10 are 162,620; 149,496; and 149,351. Three mechanical/validity tests passed.

Matched-site 3Di/amino-acid/geometry benchmarking reproduced all 858 whole-marker and 738 domain geometric baselines. Across three confidence regimes it accepted 4,743 scope/regime rows and explicitly excluded 45 for insufficient coverage. The strongest regime retains 856 whole-marker and 717 domain comparisons. Geometry is recomputed on exactly the same sites as the state mismatches. The figure was visually inspected and its RMSD/distance axes start at zero. The PH_SPT16 artifact-control case illustrates that confident local features coexist with uncertain distant-domain placement. These are dependent descriptive state fractions, not substitution-corrected branch lengths or physical conversions. Branch-estimate benchmarking, model suitability, source circularity controls and all eight evolutionary analyses remain incomplete. See docs/structural-alphabet-benchmark.md.

## Paired structural phylogenetic checkpoint

Prepared and audited all 65,750 taxon-marker cells from the 526-taxon/125-marker design. The current audited structural snapshot supplies 52 eligible paired AA/3Di alignments. Retrieved the authors' two empirical structural matrices from Edmond version 3.0 with published MD5 and local SHA256 verification. Completed 208 IQ-TREE fits and independently audited 416 matched unrooted branches on sequence-derived marker topologies. Seven focused tests passed. Warnings remain explicit: 169 fits emitted warnings and 85 contain a branch ≤10⁻⁵. A source audit explains small-tree memory allocation messages without dismissing rare-state/near-zero concerns. Full sampling, species-tree support, uncertainty, model adequacy and evolutionary objectives 1–8 remain incomplete; current branch estimates are a checkpoint, not structural-acceleration findings. See paired-structural-phylogenetics.md and metadata/paired_* for execution evidence.

## Paired uncertainty checkpoint and inventory completion

Completed 31,200 paired resampling attempts for all 52 currently covered markers under block lengths 1/10/30. The full audit verified 62,396 fits; two unestimable draws were retained without replacement. Produced 2,496 conditional interval summaries and an inspected reproducible figure. Five new focused tests passed. An independent native-feature overlap audit quantifies substantial dependencies beyond local blocks; these results do not establish calibrated acceleration tests or finish evolutionary objectives 1–8.

UniProt matching reached terminal records for all 526 working taxa: 525 completed queries (including zero matches), and Chromosphaera OFS5426494 remains taxonomy-ID unresolved. The prior 183-taxon retrieval snapshot was gracefully stopped after 16,709 additional accession attempts, draining active downloads and exiting successfully; a refreshed full-inventory queue was launched with the same two-worker retrieval policy. Model validation and broader structural coverage remain ongoing.

The refreshed queue validated 23,633 cached accessions and now prioritizes 19,209 remaining marker accessions across 196 taxa before the other atlas candidates (2,990,876 accessions queued in total). This is candidate retrieval, not a count of final usable models.

## Tree-path benchmark and compressed-model retrieval fix

Completed the leaf-path/direct-geometry benchmark for all 52 currently covered marker trees: 722 accepted taxon pairs, eight coverage exclusions, jointly propagated path intervals and an inspected four-panel figure. Two path/covariance tests passed. The automatically selected largest-RMSD pair underwent an additional Pfam-region/PAE diagnostic, demonstrating a strong competing explanation involving uncertain relative placement. These are comparative benchmarks and artifact checks; the full evolutionary objectives remain incomplete.

The full-inventory retrieval queue encountered a gzip-encoded contributed model for B9W6V0 despite a .cif URL. After a clean stop at 1,071 attempts, added gzip decoding with response provenance and an exact literal-sequence fallback for CIFs lacking the canonical sequence field. Three focused tests passed, and the real failed accession now verifies both the original GDM model and the contributed ATBC ColabFold model against independently nominated complete sequences and CA positions. No modification codes are guessed. The error record is retained and a successful receipt appended. Retrieval resumed with the same two-worker policy. These alternative pipeline models are potential source-sensitivity controls, not experimental validation.

## Explicit prediction sources and expanded GDM mapping

Started a new marker mapping across the full design using an exact provider/tool policy (GDM, AlphaFold Monomer v2.0 pipeline), a frozen inventory, and per-marker source-availability audit. Confidence scores no longer select between different pipelines. The producer's terminal receipt remains required before claiming expanded coverage.

Completed every available cross-pipeline exact-sequence comparison in that frozen inventory: one 173-residue sequence, B9W6V0, with GDM and contributed ATBC ColabFold models. Both sequences, complete CA coordinates and PAE matrices validate. At shared pLDDT≥70, 171 residues give 1.086 Å RMSD and 0.193 Å PAE-filtered local distance change. This is a limited control, not a general systematic-source estimate. Added HTTP gzip support to PAE retrieval and verified the real compressed contributed PAE file. Source-selection test passed. Expanded mapping, broad pipeline/experimental controls and the full evolutionary objectives remain incomplete.

The expanded mapping was revalidated live as PID 437405 with active CPU use and a growing residue-map output. A prelaunch PAE plan identifies 13,153 expected distinct source-eligible marker models and 412 candidate cached PAE receipts. Its conservative uncompressed-JSON estimate is 56.8 GB, with two retrieval workers and cache readback required. These are planning counts; final mapped coverage and PAE acquisition remain pending the producer's successful terminal receipt. See metadata/expanded_pae_resource_plan.json.

## Resolved final outgroup inventory query

Verified Chromosphaera perkinsii (OFS5426494) as NCBI taxon 1932427, database name Ichthyosporea sp. XGB-2017a, through an explicit NCBI alternate-name record plus the original genome BioSample SAMN06200032 in PRJNA360047. Archived ENA, NCBI and article evidence with checksums. A versioned, evidence-checked query override preserves the historical manifest and both scientific-name representations. The override integrity/conflict test passed. Its UniProt query completed with zero records; all 526 taxa now have completed query records, not universally available structures.

Prepared the full-inventory prediction snapshot: all 58,883 unique marker sequences retained, with 14,516 verified-reuse receipts, 16,106 pending nominated retrievals and 28,261 no-candidate sequences. There are no remaining inventory-pending sequences. The audited delta against the original candidate queue contains 7,799 additional sequences (4,252 short/canonical), while 18 earlier candidates now have reuse evidence. This is input preparation; the running GPU job retains its frozen inputs, and new predictions have not been submitted from this delta. See chromosphaera-taxonomy-resolution.md for sources and commands.

## Expanded catalog, confidence prefetch and lineage coverage

Validated 13,153 distinct source-specific models with 13,654 marker links across 322 taxa from the mapper's frozen inventory. Confidence-matrix prefetch is running independently with two workers after catalog validation and the existing resource estimate; final mapping identity agreement and a mapping-bound PAE receipt remain required downstream. The residue mapper was directly revalidated live as PID 437405; an incomplete output is not accepted as completion.

Completed and independently hash-checked the full 526-taxon/125-marker coverage audit. Four exclusive states contain 13,654 catalog-model links, 17,764 reuse candidates outside the catalog, 28,422 recovered marker links without a query candidate, and 5,910 unrecovered marker slots. Linked taxa comprise 304 fungal entries and 18 outgroups; 204 taxa lack catalog models. Published the per-taxon and per-lineage tables and visually inspected the figure. Separately frozen acquisition states and model-source differences remain explicit, avoiding interpretation of acquisition-order bias as biological absence. Expanded mapping, completed confidence retrieval and evolutionary inference remain pending.

## Expanded mapping completion and gene-copy annotation audit

The expanded mapping reached a successful terminal receipt: 4,895,692 residue links, 13,654 marker proteins, 13,153 models, 322 taxa, and 124 markers with at least four model-linked taxa before confidence filtering. All artifact hashes and exact model/marker identities agree with the prefetch catalog. Started four-thread native structural-alphabet extraction across 6,910,765 full-model residues after recording its resource envelope. PAE retrieval remains running; final confidence-bound native auditing is pending.

Completed the full 526-taxon raw-BUSCO/gene annotation audit. Of 2,305 duplicated calls, 312 map to one annotated gene's multiple products, 1,989 to multiple annotated genes, and four have unresolved mappings. Four original complete-single-copy hits are excluded by representative selection and remain explicit. All per-taxon raw statuses agree with archived QC; four focused tests passed. F1243177 has 119 multi-gene duplicated markers, so its extreme duplication score is not explained by annotated isoforms. Assembly redundancy, biological copy number, contamination and reconciliation still require assessment; production marker calls were not altered.

Expanded native extraction then finished successfully. All 13,153 model exports passed independent identity and finite-descriptor readback: 6,910,765 residues, 69,107,650 feature values and 26,306 zero-offset rows retained for explicit validity treatment. Output occupies approximately 866 MB. Native coordinate/partner reconstruction, completed PAE acquisition and confidence-qualified phylogenetic inputs remain pending.

## Staged native validation and expanded coordinate execution

Separated coordinate/partner reconstruction from PAE qualification without changing the final feature definitions. Coordinate-only outputs explicitly omit PAE arrays and joint-confidence counts. The confidence stage requires identical mapping provenance and complete model/version sets, verifies per-model encoding/PAE hashes, and retains directional confidence and invalid-site missingness.

Replayed the completed 422-model reference through both stages and independently compared every model summary and all 2,954 arrays: exact agreement, including NaN placement. Three geometry and two directional-context tests passed. Launched the coordinate-only audit across all 13,153 expanded models with one numerical worker after recording the resource envelope. This audit and confidence retrieval are running; expanded confidence-qualified inputs and evolutionary results remain pending. See prediction-source-controls.md for commands and the distinction between completed regression and ongoing production execution.

## Duplicated-marker sequence and genomic-location audit

Compared all 3,189 protein pairs within 2,305 duplicated BUSCO calls, verifying original proteome identities and retaining the full 526-taxon summary. Found 238 identical full-protein pairs assigned to distinct annotated genes, across 229 marker calls in 75 taxa. All copies remain separate. Resolved all pair locations from checksum-matched annotations: 208 on different assembly sequences and 30 in disjoint intervals on the same sequence, with no overlapping or unresolved intervals. One Creolimax pair explicitly uses CDS bounding spans rather than gene-feature intervals. Two source-format tests and independent hash/count/identity checks passed.

These are quality-control candidates for assembly/flanking-sequence and gene-tree review, not confirmed duplication or redundancy events. Expanded coordinate validation, confidence acquisition, local predictions, marker/species trees, orthology and full-proteome domain searches remain active. The full evolutionary objectives remain incomplete.

## Full MAFFT tree launch and taxon-coverage sensitivity

Launched an independent guide search on the complete 526-taxon, 63,750-site MAFFT matrix after validating all source artifacts and identical taxon membership against the profile matrix. Same model/seed as the profile guide, 16 threads and 32 GB memory; prelaunch resource estimate recorded. IQ-TREE entered likelihood optimization and estimated 20,661 MB RAM. All 526 sequences fail its nominal composition test, so richer-model/compositional analyses remain necessary. Neither guide is a supported final species tree.

Completed direct matrix coverage verification for all taxa and recorded 10/30/50/70-percent joint-method sensitivity memberships. These retain 522/507/495/472 taxa, respectively. A 50-percent cutoff would remove all sampled Microsporidia and Olpidiomycota; 70 percent also removes Sanchytriomycota. Production inputs retain all 526 taxa. F1243177, F1507870, Pirum and Abeoforma fall below 10 percent in both matrices, with copy/marker causes explicit. Filtered-tree inference and topology sensitivity remain pending; no placement is declared reliable solely from coverage.

## Twelve completed marker trees verified

Completed a second immutable marker-tree support audit: 12/125 trees, 5,742 internal branches, 5,630 reported support values and 112 unreported values. Exact input/tree taxa, source hashes, branch validity and output hashes pass. The remaining 113 markers are explicitly pending; no full-dataset discordance or species-tree support result is claimed. Both full-matrix guide searches and expanded structural validation remain running.

## Follow-on prediction inputs and length-aware resources

Prepared a runner-compatible queue of all 7,799 additional missing-model sequences, retaining 7,895 links across 156 taxa and excluding every original candidate. Checked current verified AFDB models and existing local predictions; no overlap required reuse at this snapshot. There are 4,252 short/canonical candidates and 3,547 explicitly deferred by length/alphabet. Input hashes and disjointness passed independent readback.

A projection from 3,640 frozen prediction receipts stratified by sequence length estimates 8.23 inference GPU-hours for the short follow-on queue, with a 16.77-hour scheduling allowance. This is preparation, not an additional GPU launch. The active ESMFold run remains on its original queue; fresh reuse and GPU checks precede the next execution. Expanded coordinate audit, confidence retrieval and phylogenetic/orthology/domain jobs continue.

## Assembly-matched coding sequences started

Started checksum-validated CDS acquisition for all 519 NCBI-backed taxa after recording a two-worker resource envelope. Every URL matches its protein assembly directory/prefix. Early downloads pass publisher MD5, local SHA256, gzip, identifier and DNA-alphabet checks; non-triplet, ambiguous and unlinked records remain explicit. Three focused tests passed. The producer is running and the completed 519-taxon receipt is not yet available.

Separately verified the CDS source inventory for all seven external taxa. Four published CDS FASTAs are present and checksum-valid, two sanchytrid translated/extracted subsets remain available with their original exception records, and Creolimax extraction is pending. This advances selection-analysis inputs; translation/codon validation and selection tests remain incomplete. See coding-sequence-workflow.md and metadata/external_cds_source_inventory.json.

## Creolimax coding-sequence extraction

Completed the genome/GTF/protein-matched extraction audit for all 8,694 Creolimax proteins. Exact translation passes for 8,558 CDSs; 136 remain exceptions (71 translation mismatches, 60 initial-phase flags, four non-triplet lengths, one ambiguous transcript mapping). Four strand/splice/phase/bounds tests passed, and independent output readback reproduced every accepted translation and the complete status partition. No exception was silently repaired.

Updated the external CDS source inventory: all seven external taxa now have available published CDSs or verified extracted subsets; published-source translation checks and extraction exceptions remain explicit. Full NCBI CDS acquisition continues, as do the structural and phylogenetic jobs. Selection tests remain pending.

## Expanded coordinate audit completed

The 13,153-model native coordinate/partner audit reached successful completion: 6,910,765 residues, 6,884,459 valid states and 26,306 invalid terminal states. Six-residue pLDDT70 retains 4,719,638 valid states. All 13,153 encoding hashes, model identity counts, summary totals and script provenance passed independent readback. Encodings occupy approximately 95 MB. PAE retrieval is still running; no joint-confidence result or expanded evolutionary inference is claimed.

## Published outgroup coding boundaries resolved

Completed strict translation auditing for all 59,116 proteins in four published outgroup proteomes: 41,592 direct passes, 17,399 non-triplet flags and 125 translation mismatches. Genome/GFF reconstruction with explicit initial-phase and terminal-remainder handling reproduces all 59,116 proteins. The 12,431 apparent genomic mismatches from the initial audit are fully explained by published CDS exports that already removed the annotated initial phase. Both original and revised audits remain preserved; no arbitrary frame search was used.

Seven focused tests pass, and independent readback verifies every output hash and accepted translation. Partial boundaries remain explicit; codon-selection eligibility and selection inference are pending. See coding-sequence-workflow.md for commands, source conventions and receipts.

## All NCBI coding-sequence downloads completed

The full 519-taxon CDS acquisition finished with zero pending/error taxa. Files occupy 2.77 GB compressed. Publisher MD5, local SHA256, gzip integrity, unique identifiers and DNA alphabet passed during acquisition; partial, ambiguous and protein-unlinked records remain counted. Independent receipt readback checks every file/hash and exact assembly-matched manifest coverage. All seven external-source taxa have separately documented CDS sources or verified derived subsets. NCBI translation validation and codon/selection feasibility remain pending.

## Full NCBI translation audit launched

Started the strict translation audit across all 519 NCBI-backed taxa and 5,847,336 CDS records, using checksum-verified GFF and protein inputs. Translation codes are annotation-derived or explicitly marked as table-1 assumptions. Duplicated mappings, partial lengths, source exceptions and mismatches remain visible. Four focused tests passed. Initial taxa completed, confirming the producer is executing; no whole-dataset translation result is yet claimed. Resource plan and restart checks are recorded in the coding-sequence workflow.

## Marker-to-CDS identity indexing started

Started source indexing for all 59,840 marker/taxon links across the full 526-taxon dataset. The index checks exact protein sequences and preserves gene-representative decisions, source-boundary conventions, missing CDSs and ambiguous duplicate records. It exports uniquely associated source DNA for later translation-qualified codon alignment. Two identity/ambiguity tests pass; initial taxa completed. This indexing run and the independent full NCBI translation audit are both active; no full codon-alignment result is claimed.

## Full marker CDS source index completed and verified

All 59,840 marker/taxon links across 526 taxa have a unique CDS source; no missing or multiply linked source records remain in this marker set. Independent readback verified the complete identity universe, original marker provenance, all exported nucleotide hashes/lengths and manifest coverage. Exported DNA totals 100,455,973 bases.

The index preserves four alternative products and 313 unresolved/provisional gene mappings, recorded in a 317-row exception table. These are identity/annotation flags, not missing sequences. Full NCBI translation auditing continues; no codon-selection eligibility or selection inference is claimed. Verification script, receipts and coverage/exception tables are versioned.

## Marker code provenance and partial-boundary audit started

Verified NCBI's documented default-code and GFF boundary conventions against its primary documentation. Started a full 519-taxon marker annotation audit to record CDS/region/default code provenance, strand-aware terminal and internal partialness, phase continuity, overlaps and strict source-CDS translation. The original full-proteome producer remains unchanged and running. Five focused tests passed; initial taxa completed. Annotation/translation matches remain separate from genome reconstruction and codon-selection eligibility.

## Full marker boundary audit completed; codon projection started

The completed 519-taxon audit covers 59,269 marker sequences: 58,967 exact translations, 280 non-triplet lengths and 22 mismatches. Ordered feature checks pass for all records, with 919 partial 5-prime ends, 965 partial 3-prime ends and 201 nonzero initial phases retained as separate flags. No internal partial boundaries or code conflicts were found. Code provenance is documented default for 58,562 markers and explicit CDS metadata for 707. Artifact hashes, full identity coverage and recounted summaries passed independent readback.

Started codon projection for all 125 markers using the unchanged full MAFFT protein alignments and 50-percent occupancy columns. Five projection tests pass. Excluded records remain explicit; accepted partial and gene/isoform flags remain attached. Independent output translation readback is prepared and will run after projection completes. These alignments are not yet qualified for selection models.

## Full codon alignment projection completed

All 125 marker projections finished: 59,334 accepted sequence rows, 506 explicitly excluded, and 63,750 unchanged protein-mask columns represented as codons. Accepted taxa range from 418 to 502 per marker. There are 123 mixed-code markers, requiring code-compatible modeling subsets rather than one indiscriminate genetic code. Producer receipts and compact review rows are versioned. Independent per-codon translation readback is running; its success is not yet claimed. Selection feasibility, model fitting and saturation assessment remain pending.

## Codon readback completed

Independent verification passed for all 125 codon alignments and all 59,334 accepted sequence rows. Every one of 28,328,533 non-gap codons translates under its recorded code to the exact original masked amino-acid column. Full source identity coverage, artifact hashes, triplet gaps, absence of aligned stops and marker totals pass. The 506 exclusions remain explicit. Receipt and marker summaries are versioned.

Started the full marker coverage screen for groups sharing fungal genus labels and translation codes, comparing all aligned entries with exclusion of recorded annotation/gene/label flags. Three coverage tests passed. Group labels are not assumed monophyletic or independent ecological transitions; this screen does not estimate selection or establish divergence suitability.

## Genus/code codon coverage screen completed

Screened all 125 markers across 16 fungal genus labels containing 113 manifest entries. Of 4,168 marker/code/policy rows, 1,826 groups pass the inclusive coverage threshold and 1,712 pass after excluding recorded annotation/gene/taxon-label flags. Output hashes, row identity, retained membership, stricter-policy subsets and threshold totals pass readback checks.

Serendipita loses all passing groups under the stricter policy because three of its five entries are unnamed species records; only two named entries remain. This is a documented taxon-identity constraint, not a finding of missing sequence or ecological nonreplication. Fifteen labels retain passing groups under the stricter policy. Genus monophyly, divergence/saturation, selection-model adequacy and independent ecological transitions still require assessment.

## Full NCBI translation audit and readback completed

All 519 taxa completed: 5,847,336 CDS records and 5,843,347 normalized proteins, with no protein lacking a source CDS association. Exact translation passes for 5,752,457 records. The remaining exclusive classifications comprise 85,683 non-triplet lengths, 2,491 mismatches, 2,715 annotation exceptions, 3,988 missing/ambiguous protein identifiers and two CDS records sharing a protein.

Readback checked all output hashes and recounted every CDS/status/code total and taxon identity. All 59,269 NCBI marker classifications reconcile with the separate boundary audit: 59,252 identical statuses and 17 differences explained by the full audit's annotation-exception gate. No unexplained code/status disagreement remains in the marker subset. Full nonmarker coordinates, region-code handling, family codon suitability and selection inference remain pending. Receipts and taxon summaries are versioned; large per-CDS tables remain outside Git.

## Observed codon divergence and review figure completed

Measured shared-site observed differences for 43,235 distinct marker/code/taxon pairs across all coverage-passing candidate groups. The two inclusion policies produce 84,796 rows; repeated metrics agree exactly. Source/output hashes, numerator partitions, fractions and overlap thresholds pass consistency checks; four focused tests passed.

Rendered and inspected a figure of all 1,712 stricter-policy marker/group medians. The highest median amino-acid difference is 0.7371 for Aspergillus marker 4986044at2759 despite at least 232 shared codons per pair. Twenty leading cases are flagged for alignment/orthology review, not biological acceleration. These observed fractions do not estimate dS, dN/dS, saturation or branch-specific rates. Independent group realignment and profile-alignment comparisons remain next checks before selection modeling.

## Outlier realignment and mixed-domain evidence

Completed all 20 targeted group realignments and 163 pair correspondence comparisons, reproducing every original divergence baseline independently. Three coordinate/codon-mask tests passed. High median differences persist, while residue correspondences and selected columns show sensitivity to alignment context.

The leading Aspergillus marker has two domain-hit patterns: six longer BRF1-hit proteins and four shorter proteins without that hit, both carrying TFIIB hits. The global median difference is 0.75 across patterns but 0.1336 and 0.0460 within them. All ten source identities and sequence hashes, domain source hashes and output artifacts were verified. This is exploratory evidence for mixed protein types/hidden paralogy as a competing explanation, not proof of duplication history or acceleration. Gene-family reconstruction/reconciliation is now necessary before biological interpretation of this case; baseline marker inputs remain preserved.

## Broad candidate-family search launched

Verified the protected core orthogroup table and traced the only focal Aspergillus core taxon to OG0001567. Nine focal taxa are outside the reference core, so core membership alone cannot resolve the mixed-domain case.

Started a focused three-profile Pfam search across all 5,654,720 additional representative sequences from the full 526-taxon design, with existing marker hits available for later reuse. Profile/library/input/marker-source and HMMER-binary hashes pass. The one-process/two-worker search is confirmed running, with resource estimates and configuration versioned. Full homolog recovery, family phylogeny and reconciliation remain pending; original full-Pfam and orthology jobs continue independently.

## Focused full-proteome family search completed

All three selected Pfam profiles finished searching the 5,654,720 additional unique sequences in 64.54 seconds, returning 2,109 domain-hit rows. All output and producer hashes, row counts and profile identities passed readback. This completes focused acquisition for the additional-sequence partition; combining reusable marker hits and expanding sequence-to-gene/taxon mappings is still required before family inference.

## Candidate homolog identities restored across the full design

Combined focused and marker domain hits to recover 1,485 representative protein entries (1,434 unique sequences) across 524 taxa, preserving every taxon/gene identity. All exported sequence hashes/lengths, source artifacts and gene mappings passed validation. Two taxa have no candidate in this profile screen; no biological loss is inferred. The candidate union remains domain-defined, with varied domain-copy counts and eight unresolved/provisional gene mappings.

All ten focal Aspergillus genomes contain both domain-hit types on distinct annotated genes, but the original marker selects one type in six taxa and the other in four. The 20-copy identity table is versioned. This identifies copy-sampling inconsistency as a concrete explanation to test before interpreting divergence as acceleration. Domain-matched family alignment, supported gene trees and reconciliation remain next steps.

## Domain-matched family alignment completed; supported tree started

Aligned 1,040 proteins from 512 taxa using two ordered, nonoverlapping TFIIB domains with at least 70-percent profile coverage each. The resulting 184-column matrix retains 445 excluded candidates with explicit hit-count/coverage/overlap reasons. Three focused tests passed. Independent readback verified every alignment/mapping cell and 175,699 source-matched non-gap residues.

Started the full candidate-family IQ-TREE search with restricted model selection, SH-aLRT and NNI-refined ultrafast-bootstrap support (1,000 replicates each), four threads and 8 GB memory. The process entered likelihood optimization successfully. Full tree completion, repeat-specific sensitivity, duplication rooting and species-tree reconciliation remain pending; this is not a final species tree or selection result.

## Repeat-specific family sensitivity launched

Prepared two 92-column matrices with exactly the same 1,040 genes as the combined domain alignment. Independent split readback reproduces every original sequence; per-repeat coverage, identical-sequence counts and site-variation summaries are recorded. One second-repeat column has no canonical residues and remains explicit.

Both repeat-specific supported tree searches entered likelihood/model optimization successfully, with two threads and 4 GB memory per job. Each requests 1,000 SH-aLRT and 1,000 NNI-refined ultrafast-bootstrap replicates. The combined-domain tree remains running. Topology/support concordance and duplication/species-tree reconciliation remain pending; no repeat-history conclusion is claimed.

## Expanded PAE acquisition completed; confidence integration running

All 13,153 catalog models now have verified PAE downloads, with zero failures: 1,779,429,556 compressed bytes and 12,518,772,192 original JSON bytes. The completed manifest hash passes readback. The catalog's models correspond to the expanded final mapping with 13,654 marker links across 322 taxa; these counts describe structural availability, not confidence-qualified phylogenetic coverage.

The dependency-gated controller observed termination of the original prefetch process, checked its full completion receipt and exact catalog/mapping agreement, and started cache revalidation against the final mapping receipt. It will then qualify the completed coordinate encodings using six-residue directional PAE maxima and independently recount per-model joint pLDDT/PAE filters. Configuration, producer hashes and resource estimates are versioned; both directional-context tests pass. Final confidence qualification and expanded paired evolutionary analyses remain pending. The three TFIIB family/repeat tree searches continue.

## Expanded paired evolutionary input preparation queued

Started a live dependency-gated job for the expanded paired amino-acid/3Di alignments. It waits for the actual confidence-controller process to finish and requires its verified completion receipt before preparing inputs. All 526 taxa and 125 markers are screened using the same 49,027-column profile matrix and the same joint confidence/coverage rules as the earlier snapshot. Five paired-input/model-validation tests pass; pinned configuration and a one-worker, 8 GB planning estimate are versioned.

After preparation, the job will compare per-marker eligible-taxon counts and retained columns against the checksum-identified earlier snapshot. Prediction availability and selected models can change, so coverage differences will not be attributed solely to sample size. No expanded paired alignment completion, likelihood fit or acceleration result is claimed yet.

## Expanded branch-fitting preparation

Updated the idle paired-fitting runner to derive dataset dimensions from the
completed input summary and optionally retain 1,000-replicate sequence-topology
SH-aLRT/NNI-refined ultrafast-bootstrap support. All bootstrap tree identities
and counts are checked before accepting a fit. Three new support-control tests
and two existing branch-correspondence tests pass. Structural branch estimates
remain conditional on the sequence topology; branch intervals and topology
sensitivity are separate analyses. Expanded fitting is not launched while
confidence integration and paired-input preparation remain incomplete.

## Dedicated predictor-control inputs completed

Prepared 266 exact-sequence AlphaFold references for later ESMFold prediction,
covering 126 taxa, 23 lineage/role groups and 78,593 residues. Selection balances
available length/confidence strata within manifest groups; all source model and
taxon/marker identities remain explicit. Output and coordinate hashes and
sequence identities pass independent readback. Resource estimates use current
production timings (0.467 GPU hours of median-bin inference, 1.10 hours with
allowance). Prediction and geometry comparison are pending availability of the
currently occupied authorized GPU; the missing-structure production continues.

## Paired lineage coverage limitation quantified

Completed the full 65,750-cell taxon/marker coverage summary for the earlier
paired snapshot. Independent checks against all 52 actual paired FASTA sets
show that the 105 retained taxa are entirely Ascomycota, with only two retaining
at least ten markers. All 27 manifest lineage/role groups remain in the output,
including zero-coverage groups. Prior branch results cannot support
cross-phylum structural-rate claims. Expanded lineage coverage is pending its
completed paired inputs.

All 13,153 PAE files now also pass cache revalidation bound to the final expanded
residue mapping, with zero failures. The exact manifest and mapping receipt
hashes pass readback. Confidence qualification is running; its joint-filter
counts are not yet claimed.

## Dedicated predictor controls launched on newly idle GPU

The previously occupied GPU 0 is now independently verified idle. Started all
266 controls there using the existing pinned ESMFold environment and settings;
GPU 1 continues the missing-model queue. Both live processes were observed on
their separate cards, and the first four control predictions completed. The
new configuration differs only in input receipt and GPU UUID. The launch
resource record supersedes the earlier waiting assumption without adding paid
resources. Extended independent artifact auditing to accept the reference-link
format, with three identity-validation tests passing. Full control completion,
readback and predictor geometry comparisons remain pending.

## Predictor geometry comparison implemented

Prepared exact-sequence AlphaFold–ESMFold comparison with complete-prediction
and independent-readback gates, whole-protein CA RMSD, local distance changes,
joint pLDDT thresholds and bidirectional PAE sensitivities. Four focused tests
pass. Resource estimates and the reproducible command are versioned. The
control GPU job remains live; no completed empirical comparison is claimed.

## Genus-group phylogenetic checks executed on completed markers

Audited 17 completed marker trees (108 remain pending) and 8,031 internal
splits. Across the 16 codon-screen genus labels, 250 of 272 marker/group rows
have sufficient taxon coverage; 20 contain an incompatible split with SH-aLRT
support at least 80. Naganishia shows conflicts in five of ten assessable
markers. Missing group members and all results remain explicit. Two tests and
output/hash/full-grid checks pass. This partial gene-tree screen prioritizes
orthology/alignment review and does not establish species monophyly, selection
or ecological replication.

## Expanded structural confidence integration completed

All 13,153 models passed the complete staged pipeline and per-model readback:
6,910,765 residues, 6,884,459 valid native states and 4,714,151 states meeting
joint six-residue pLDDT ≥70 / maximum directional context PAE ≤10. Both source
linkage and the final summary hash pass. Whole-model state counts do not equal
aligned or phylogenetically eligible sites. Completed receipts are versioned.

The dependent paired-alignment producer has started on the full 526-taxon,
125-marker matrix using these qualified encodings. Lineage coverage and
expanded likelihood fits will follow its completed outputs. The overall
structural atlas and evolutionary analyses remain incomplete.

## Expanded paired alignments and lineage coverage completed

All 526 taxa and 125 markers were screened. The expanded snapshot produces
paired alignments for 124 markers (5–135 taxa, 73–1,184 columns), with 322 taxa
in the union: 304 fungal entries and 18 outgroups. All emitted AA/3Di identity
sets, dimensions and masks pass readback, and actual file memberships agree
with every per-taxon coverage count. Twenty-three lineage/role groups are
represented; Aphelidiomycota, Calcarisporiellomycota, Sanchytriomycota and
Corallochytrea remain uncovered. Source entries remain in the full design.

Launched the expanded 496-fit batch with four one-thread workers, 2 GB per fit,
sequence-tree SH-aLRT and NNI-refined ultrafast-bootstrap support, and three
structural models conditional on each AA topology. The first marker is in
likelihood optimization. This begins expanded branch estimation; completion,
uncertainty/model adequacy and biological interpretation remain pending.

## Expanded lineage coverage figure completed

Rendered and inspected the paired-coverage figure across all 27 lineage/role
groups. Depth remains skewed: 204 entries have zero usable markers, 200 have
1–9, 17 have 10–49 and 105 have at least 50; all of the last group are
Ascomycota. This limits interpretation of broad group presence as sufficient
cross-lineage information. The figure, reproducible script and source/artifact
hashes are versioned. Ongoing acquisition and missing-model prediction remain
necessary to improve depth outside the best-covered lineage.

## Active ESMFold production snapshot readback started

Started independent validation of a frozen list of 5,121 per-model receipts
from the ongoing missing-structure production. The new explicit snapshot mode
preserves partial status and does not mistake an old chunk receipt for current
queue completion. Two integration tests and three link-identity tests pass.
The audit process is live. Completed artifact validation, lineage coverage and
ESMFold residue/structural mapping remain next steps; original production and
the separate 266-sequence predictor controls continue on their GPUs.

## Audited ESMFold snapshot completed and conversion started

All 5,121 frozen prediction models passed readback, linking 201 taxa and 86
markers. The subset includes both sampled Aphelidiomycota, but does not yet
establish usable structural coverage there. Started exact-field PDB-to-mmCIF
conversion with canonical polymer sequences and local prediction provenance;
over 500 models have passed actual conversion roundtrips. Three conversion
tests and three mapping-source/identity tests pass. Mapping now supports a
separate local inventory without manufacturing UniProt identities. Full
conversion, residue mapping and native confidence qualification remain pending.

## Local prediction PAE integration prepared

Implemented mapping-bound local PAE export with original prediction/NPZ
provenance and exact matrix roundtrip checks. Three focused tests pass. The
exporter explicitly records local derivation, rather than remote AFDB retrieval,
and keeps asymmetric directional values intact. Launch awaits completed local
mmCIF conversion and residue mapping. Conversion and the separate control
prediction run remain active.

## Predictor controls completed; follow-on production launched

All 266 controls completed and passed independent artifact audit. The full
798-row AlphaFold–ESMFold geometry comparison is complete: 209 proteins pass
joint pLDDT ≥70 coverage, with median RMSD 1.114 Å and local distance change
0.236 Å. Same-protein cohort summaries separate cohort changes from residue
selection effects. Results describe predictor agreement, not experimental
accuracy or biological acceleration.

Launched the 4,252 eligible follow-on missing-model predictions on the newly
free GPU 0 after confirming disjoint inputs and no existing AF cache matches.
The first prediction completed. Original prediction, local mmCIF conversion,
phylogenetic and expanded paired branch-fitting jobs continue.

## ESMFold conversion completed and marker mapping launched

All 5,121 converted models passed inventory and per-file checksum readback,
following exact-field conversion roundtrips. Started the local ESMFold-only
mapping against the full profile-alignment/matrix source design. Mapping,
local PAE binding and native feature qualification must finish before counting
new confidence-qualified coverage. The 124-marker GDM branch batch, full
phylogenetic searches and both missing-model GPU queues continue.

## Predictor agreement figure and inspection order completed

Rendered and inspected the same-209-control comparison figure. It distinguishes
whole-protein versus confidence-filtered RMSD from confident local distance
changes, with excluded controls and residue retention explicit. All qualifying
controls are retained in the global-RMSD inspection order. Large global/local
differences motivate domain/orientation review but do not establish accuracy,
biological acceleration or a causal explanation. Figure, script, table and
provenance are versioned; local ESMFold marker mapping remains running.

## Local marker mapping verified; PAE and native audits advancing

Completed the ESMFold-only mapping: 5,121 models, 5,161 links, 201 taxa and
900,303 mapped residues. Every mapped sequence/matrix character and confidence
was checked against the original local predictions, and link identities exactly
match the independent artifact audit. Started provenance-bound PAE export.
Native Foldseek extraction finished, and coordinate feature reconstruction is
now running. Joint-confidence eligibility and improved paired lineage coverage
remain pending; no new branch or accuracy conclusion follows from extraction.

### Local confidence-qualified paired coverage completed

Completed local PAE export and coordinate qualification for all 5,121 frozen
models; independent encoding readback confirms 842,713 joint-confidence states.
Prepared and independently checked 72 local paired marker alignments with
4,669 usable taxon-marker combinations across 199 taxa. Comparing the separate
prediction-source inputs increases represented entries from 322 to 458
(440 fungal entries, 18 outgroups), including 136 newly represented entries.
Amoeboaphelidium occidentale contributes 26 usable markers, adding paired
coverage for Aphelidiomycota. Coverage remains uneven and source-specific;
these counts are not completed branch tests or a fully sampled atlas.

Both ESMFold production queues and the expanded supported paired fits remain
active. Species/gene trees, orthology assignment and whole-proteome domain
searches also continue. The overall project goal remains active and incomplete.

### Local supported branch fits and predictor-alphabet control

Started 288 fits across 72 local marker alignments (4–158 taxa; 73–422 sites):
AA topology search with 1,000 SH-aLRT and 1,000 UFBoot replicates/NNI refinement,
then three structural models on each fixed AA topology. Four single-thread
workers, 2 GB per fit, existing host only. Resource estimate and launch
configuration are versioned; outputs remain incomplete.

Completed same-sequence structural-alphabet comparison for all 266 predictor
controls. At joint six-residue pLDDT70/PAE10, 187 pass coverage and show median
14.33% state disagreement and 18.18% spatial-partner changes. These are
predictor effects on identical sequences, not evolutionary change. Unit tests
and independent comparison-table/count readback passed. This result makes
prediction-source sensitivity a material requirement for branch conclusions.
The overall goal remains active; supported species trees, reconciliation,
uncertainty and the full evolutionary analyses are not complete.

### Direct local geometry and isolated orthology alignment diagnosis

Started exhaustive direct CA geometry comparisons for the frozen local snapshot:
306,673 within-marker taxon pairs, three focal-confidence thresholds, at most
920,019 comparison/exclusion rows. One CPU worker with one BLAS thread and a
16 GB planning memory allowance; existing host only. This supplies a physical
geometry baseline for structural-alphabet benchmarking, not an additive branch
metric. Joint feature/PAE masking and tree-path uncertainty are separate stages.
The producer now records its script checksum in the final receipt.

The full orthology process is still live but its log reports a FAMSA exit 134
for OG0001522 and a missing downstream alignment. Its 769 unique input IDs have
nonempty protein sequences of length 79–4,623. Started isolated reruns with the
bundled FAMSA 2.2.3 and existing FAMSA 2.5.2, each using one thread and separate
outputs. Both were confirmed running. No live OrthoFinder artifacts were changed;
no successful repair or completed orthology inference is claimed. Executable
hashes, commands, the input checksum and resource estimates are versioned in
`metadata/orthology_family_alignment_failure.json` and
`metadata/orthology_family_diagnostic_resource_plan.json`.

Updated the README opening status to reflect completed local coverage and the
same-sequence predictor-alphabet control, replacing the obsolete first-chunk
status. The preceding goal turn made progress through completed predictor
comparisons and launched supported local fits; the overall goal remains active.

The isolated OG0001522 reruns subsequently both finished with exit 0. Independent
readback verified all 769 input identities and complete ungapped protein
sequences. Bundled FAMSA yielded 20,197 alignment columns; FAMSA 2.5.2 yielded
19,540. The original crash did not reproduce and its cause remains unresolved.
Using the successful **same bundled executable** output, started isolated
production-equivalent trimming and FastTree recovery. The recovery script
checks every trimmed column against the ordered source alignment and requires
the final tree to retain every gene ID with finite nonnegative branches.
Production alignment/tree paths remain absent; integration is still pending.

OG0001522 recovery subsequently completed: 2,112 retained alignment columns and
all 769 tree tips passed validation. Verified the main job was still before
species-tree inference, with no active family writer and both production paths
absent. Atomically installed independent copies of the recovered alignment and
tree; final hashes match the isolated recovery and protected core artifacts
remain unchanged. No existing output was overwritten. The full orthology job
continues; this repair does not establish completed reconciliation.

### Full-ingroup complementary lineage QC queued

Prepared explicit selection and offline BUSCO execution for all 501 fungal
entries, preserving broad results for every taxon and retaining 25 outgroups
outside fungal-panel scoring. Acquiring six frozen OrthoDB12.2 panels and
verifying publisher archive checksums, configurations and file hashes. Five
matching phylum groups use their named panel; remaining fungal entries use the
fungal parent panel. A path-finalization defect was corrected with original
source preserved and explicit verified-archive reuse. The controller is running
acquisition/finalization and will start the full batch after successful checks.

Three summary-validation tests pass. Resource plan, controller configuration,
recovery provenance and reproducible scripts are versioned. Dataset finalization
and per-taxon QC remain incomplete; no lineage-specific completeness results
are claimed yet. See `docs/lineage-completeness-workflow.md`. The preceding goal
turn made progress by recovering the failed orthology family and launching
exhaustive local geometry; the full research objective remains active.

### Functional-site projection completed; lineage QC executing

Completed Pfam active-site projection for 10,570 hits in 8,961 marker proteins
across 173 families. Preserved 19,429 reference-pattern/site rows, including gaps,
substitutions and overlap ambiguity. The conservative candidate screen yields
5,251 distinct positions in 3,503 proteins across 94 families. All 15,704 mapped
positions passed full-protein readback, and whole-pattern/candidate statuses
were independently recomputed. Three focused tests pass. These are candidate
functional correspondences, not validated catalytic sites or evolutionary tests.
Structural and phylogenetic joins remain pending.

All six lineage-BUSCO datasets finished checksum and file verification, and the
501-entry batch has started producing successful per-taxon results. Dataset
marker counts, exact assignments and launch configuration are now versioned.
The previous goal turn advanced the executable lineage-QC pipeline; this turn
completes functional-site mapping and validates the transition to QC execution.
The overall research goal remains active and incomplete.

### Functional sites linked to structures and paired phylogenetic characters

Completed both predictor-specific joins, retaining all conserved, nonconserved,
gapped and overlap-ambiguous Pfam correspondences. Each source has the same
17,105 taxon–marker–site universe. There are 2,713 observed AlphaFold rows and
1,189 observed local rows in the actual paired alignments, spanning 399 taxa in
union. Of these, 1,258 are flagged as conserved candidates. The counts describe
coverage and repeated homologous-site observations, not independent changes.

All 3,902 observed rows passed independent native-confidence and emitted
AA/3Di-character readback. Tests passed for keeping nonconserved/gapped rows,
deduplicating reference annotations and rejecting conflicting projections.
Complete tables remain outside Git, with receipts and observed subsets versioned.
Branch localization and site-specific evolutionary tests remain outstanding.
The previous goal turn completed the initial site projection; this turn advances
its connection to the structural and phylogenetic data. Overall goal stays active.

### Species ecology expanded; prediction-source confounding measured

Expanded primary-literature curation from six to 21 species statements, retaining
both Amanita genus-projection conflicts and explicit taxonomic/isolate uncertainty.
Bound host evidence to all ten sampled Suillus entries: eight reported assignments,
one source-uncertain assignment and one unverified table mapping. No selected
isolate is marked experimentally verified or eligible for confirmatory testing.

Rebuilt evidence tables and joined both frozen predictor-specific paired inputs.
All 21 species have some eligible coverage, but the five-member Amanita group has
no marker eligible in every member within either predictor; the three-member
Cenococcum group has one such AlphaFold marker and none from ESMFold. Its
saprotrophic comparators currently have no eligible ESMFold markers, exposing
prediction-source confounding that requires matched-source coverage before tests.
Builder checks verify source hashes and exact host-panel membership. This is an
evidence/availability checkpoint, not an ecological effect or completed project.

Prepared 675 additional same-method ecological marker predictions across 18
species, disjoint from existing queues. Inputs and dispositions are frozen and
read back successfully; execution remains pending GPU availability and a resource
forecast. The exhaustive local geometry producer has now completed and awaits
independent output validation.


### Local geometry validated and ecological prediction resources estimated

Completed all 920,019 pair/threshold decisions across 306,673 local-model taxon
pairs. All eligibility and sequence-difference readbacks pass; independent SciPy
geometry matches on 223 deterministic rows. There are 677,982 accepted comparisons
and 242,037 exclusions, retained separately. Descriptive threshold summaries and
hash-pinned receipts are versioned; branch-rate benchmarking remains outstanding.

The 675-input ecological queue has a frozen forecast from 6,120 observed timings:
1.29 GPU hours median-based inference, 2.60 hours planning allowance, 24 GB VRAM
and 20 GB output headroom. Both GPUs are confirmed occupied by existing batches;
no additional device or paid resource was provisioned. The preceding turn made
progress by freezing these inputs; this turn validates completed geometry and
makes the queue's resource requirements concrete. The full goal remains active.


### Direct geometry on paired inference sites launched

Extracted direct coordinate/PAE calculations from the tree-path benchmark so the
full 72-marker local dataset can be prepared while supported tree fits continue.
All shared fields matched the prior 52-marker benchmark: 722 accepted pairs and
eight coverage exclusions. This regression verifies the refactor's behavior, not
an independent validation of its shared geometry implementation.

Launched the full local paired-site pass (up to 266,593 pairs) on one existing CPU
thread after recording memory, storage and runtime estimates. Exact paired AA/3Di
observations determine matched residues; the local metric additionally requires
both directional pairwise PAE values <=10 in each model. Input/coordinate/PAE
hashes are checked. The expanded calculation is running and has no completion
receipt yet; tree-path joins, joint uncertainty and biological inference remain
outstanding. The prior goal turn validated the broad local geometry baseline;
this turn advances comparability to the actual phylogenetic character sets.


### Ecology predictions scheduled after existing GPU work

Started the controller for all 675 frozen ecology prediction inputs, with exact
predecessor process/configuration checks, clean-chunk completion gates, device
availability checks and pinned producer/data/resource receipts. Two tests pass,
including rejection of incomplete, interrupted, OOM-deferred and mismatched
predecessor chunks. It is waiting on the confirmed live follow-on batch; no new
GPU job is claimed yet and no paid resources were requested.

The paired-site geometry and supported tree jobs remain live. The previous turn
made progress by launching comparable direct geometry; this turn closes the
execution dependency for the ecology coverage gap. The full research goal stays
active; scheduling does not complete the predictions or ecological tests.


### Hybrid evidence linked to selected assemblies; taxon sensitivities prepared

Verified accession/strain/publication linkage for VIN7 and CBS 1483 using the
frozen assembly catalogue and primary genome studies. S. pastorianus is an
additional curated hybrid missed by the name-only screen. Catalogue assembly
representation does not establish biological ploidy or nonhybrid ancestry.

Prepared four exact sequence-subset matrices: both alignment methods with the
two hybrids excluded (499 fungal entries +25 outgroups), and with those hybrids
plus 21 uncertain labels excluded (478 +25). All retained sequences passed full
readback; original data remain intact. These are unrun sensitivity inputs, not
a changed sampling objective or accepted unique-species count. The previous
turn scheduled ecology predictions; this turn adds assembly-linked evidence
and executable taxon-sensitivity inputs. The full goal remains active.


### Taxon-identity sensitivity guide trees executing

Submitted all four profile/MAFFT and hybrid/uncertain-label sensitivity matrices
to two concurrent eight-thread guide searches. Source hashes, resource envelope,
model, seed and executable identity are pinned. Expected-tip and branch-validity
checks gate completion; IQ-TREE checkpoints and verified complete-result reuse
support recovery. The existing CPU host has sufficient observed capacity.
Neither exclusion policy removes an entire role/major-lineage group, although
within-group and marker coverage remain relevant. These are running homogeneous
model guides, not supported final species trees. The preceding turn prepared
these matrices; this turn begins their actual inference. Overall goal stays active.


### Identity-policy effects on paired marker availability quantified

Reapplied the two taxon-exclusion policies to both frozen structural sources and
verified ready-family identities against actual AA FASTAs. The stricter policy
retains all 124 AlphaFold markers (320 taxa, 13,510 cells), but leaves 71 ESMFold
markers (184 taxa, 4,286 cells). Local marker 776280at2759 falls from four taxa to
three when Ceratobasidium sp. AG-Ba is removed. This is an availability effect,
not inferred structural or taxonomic change. Marker and lineage counts plus
hash receipts are versioned. The preceding turn launched taxon-sensitive guide
trees; this turn measures their structural-data coverage consequences. Existing
tree and geometry jobs remain active, and the research goal is incomplete.


### Exact TFIIB input tips annotated for copy-selection review

Joined all 1,040 tips to gene/protein, domain and BUSCO-selection evidence; all
source domain-count triplets read back successfully. The inputs have 473
BRF1-detected and 567 undetected entries, with both in 446 taxa. Selected focal
BUSCO copies include 431 detected and 33 undetected entries; 13 of the latter
taxa also have an aligned detected candidate. This extends the copy-selection
review beyond the initial genus example but does not resolve orthology or
establish domain loss. Combined/repeat trees remain live. The prior turn measured
identity-policy coverage; this turn prepares exact family annotations needed for
interpreting a mechanistic case. The full goal remains active and incomplete.


### Experimental structure candidate inventory completed

Queried the complete 13,153-accession frozen AlphaFold marker-model inventory
against RCSB experimental polymer entities. All 132 batches completed, nominating
3,592 unique entities. Saved query/response hashes, counts and the deduplicated
union passed independent readback. No atomic structures or exact-sequence
benchmarks are claimed yet; entity sequence, construct, observed-coordinate,
quality and training-overlap checks remain. The prior turn annotated TFIIB tree
inputs; this turn establishes a concrete experimental-reference candidate set
for the project's prediction/circularity controls. The full goal remains active.


### Experimental entity and entry metadata retrieval started

Launched complete retrieval for all 3,592 nominated experimental entities and
1,635 parent entries, with identity checks, atomic cached responses, hashes and
explicit resource limits. Initial entity records are arriving and their declared
identities/sequence fields are present. The full 5,227-response stage remains
running; exact sequence, construct and experimental-quality results are pending.
The preceding turn completed candidate nomination; this turn acquires the data
needed to determine actual reference eligibility. The full goal stays active.


### Experimental sequence correspondence screen executed

Implemented and tested deposited-sequence classification, then froze a partial
screen of 666 verified entity responses. Found 541 exact full canonical-sequence
matches and 66 uniquely positioned exact fragments; 59 other rows retain explicit
construct, variant/alignment or noncanonical review states. All exact-full rows
passed hash/length readback. Repeated entities are not independent proteins and
no candidate is yet declared an eligible experimental benchmark. Full metadata
retrieval continues. The prior turn started that acquisition; this turn produces
actual sequence correspondence evidence. The full research goal remains active.


### Experimental coordinate candidate acquisition completed

Downloaded all 215 entries nominated by the frozen exact-sequence screen:
324,469,519 compressed bytes and 1,297,265,050 uncompressed bytes.
Every archive passed gzip integrity and data-block identity checks; complete
entry membership, receipts and compressed hashes passed independent readback.
Two integrity tests pass. Full complexes are preserved. Atomic residue mapping,
experimental quality and independent-benchmark eligibility remain outstanding.
The prior turn established sequence correspondence; this turn acquires the
coordinates needed for actual structural benchmarking. Full metadata retrieval
continues and the overall research goal remains active.


### Experimental residue-correspondence mapping executing

Launched all 215 coordinate entries through sequence-checked entity/model/chain
mapping. Every full-sequence position remains explicit, including missing and
ambiguous CA observations; no occupancy/alternate choice is hidden. Two focused
tests pass. Initial files produce verified sequence correspondence and explicit
missing coordinates. Complete mapping, readback and actual predictor comparisons
remain outstanding. The previous turn completed archive acquisition; this turn
starts the residue mapping required for meaningful experimental benchmarks.
The full research goal remains active.


### Combined TFIIB tree completed; mixed-copy marker concern strengthened

The 1,040-tip combined-domain tree and all 1,000 bootstrap trees completed.
Independent tip/branch validation and split recount identify a 100/100-supported
split separating 555 BRF1-undetected entries from 473 detected plus 12 undetected
entries. It occurs in every saved bootstrap tree. Selected BUSCO copies fall on
both sides (437 and 27), supporting a gene-copy reconciliation requirement before
confirmatory single-ortholog structural-rate interpretation. The 48 unlabeled
internal groups contain identical aligned sequences and remain unresolved.

Two split-validation tests pass. An explicit marker caveat is versioned without
changing running/frozen fits. Repeat-tree sensitivity, rooting and reconciliation
remain pending. Experimental acquisition/mapping jobs continue. The preceding
turn launched CA mapping; this turn completes a substantive family-tree assessment.
The full research goal remains active and incomplete.


### Orthology caveat applied to current paired-marker review

Created a hash-pinned review overlay for all 250 marker/source combinations.
The flagged 116-taxon AlphaFold marker contains 108 exact protein matches on the
BRF1-enriched side, four on the other supported side and four outside eligible
family-tree coverage. Its local counterpart has only one eligible taxon. The
expanded input therefore directly exhibits the copy-mixture issue; confirmatory
single-ortholog interpretation is withheld pending reconciliation. Original fits
remain diagnostic and unflagged markers are not certified valid. The previous
turn completed the family-tree assessment; this turn connects that evidence to
the actual downstream inputs. The overall project remains active.


### Full experimental sequence screen and initial coordinate audit completed

Verified all 5,227 downloaded metadata responses and completed the 3,592-entity
sequence screen. The 2,850 full exact matches cover 80 model proteins and 1,036
PDB entries; the 215-entry coordinate snapshot remains partial coverage.

Completed and audited all 157,670 exported CA position rows in that snapshot:
129,095 unambiguous full-occupancy positions, with missing/ambiguous/modified
observations retained. Full grids, sequence hashes, embedded identities and
partitions passed; raw atom fields passed independent checks at 2,823 positions
in five entries. Experimental-quality filtering and actual prediction comparisons
remain pending. The preceding turn applied the orthology caveat; this turn
completes substantive experimental-reference acquisition and correspondence work.
The full project remains active.


### Experimental methodology review and full coordinate expansion

Completed the full reference metadata inventory and independent method,
methodology and resolution readback for all 1,036 exact-match entries. Four are
integrative models, explicitly deferred; 1,032 experimental entries remain.
Started full coordinate retrieval with verified reuse of 215 cached archives
and 817 new requests. Metadata, tests, resource estimates and provenance are
versioned; coordinate retrieval is still running.

The 72-marker local ESMFold paired-site geometry producer also finished:
259,780 accepted and 6,813 excluded taxon pairs. Output checksums passed readback;
this checkpoint does not claim a fresh independent numerical audit of all pairs.
These measurements use the same observed sites as the paired sequence/structural
alphabet inputs. Joining them to completed supported tree fits and uncertainty
estimates remains pending. Prediction, phylogeny and other background stages
continue; the overall project remains active and incomplete.


### Paired-site geometry audit completed; expanded AlphaFold run started

This continuation independently verified all 266,593 local ESMFold pair records
and recalculated geometry and directional PAE filtering for 72 deterministically
selected pairs, one per marker. All checks passed. Started the same geometry
pipeline on 124 expanded AlphaFold markers with 744,853 possible taxon pairs,
after recording resource estimates. The first marker has completed; the run
continues. Full experimental-coordinate acquisition also continues. Supported
phylogenies and downstream branch/uncertainty integration remain pending; the
full goal remains active. The preceding turn was substantive progress through
metadata review, coordinate expansion and versioned results.


### Initial experimental geometry comparisons completed

Completed 1,803 accepted prediction–experimental-chain comparisons across the
original 215-entry coordinate snapshot, with 96 explicit coverage exclusions.
Every accepted geometry result agreed with independent rotation/distance
recalculation. Summaries now report both comparison weighting and equal protein
weighting, plus fixed chain/model cohorts across confidence thresholds. The
unfiltered median is 0.900 Å per comparison versus 3.235 Å across protein medians,
showing why deposition multiplicity must be accounted for. These remain partial,
descriptive agreement results pending experimental quality/context/training
review. Full coordinate retrieval and expanded phylogenetic/structural jobs
continue. The preceding continuation completed the paired-site audit and started
expanded AlphaFold geometry; both turns made concrete progress. The full goal
remains active and incomplete.


### Full experimental archives verified; expanded residue mapping running

Finished all 1,032 coordinate archives (1.95 GB compressed), with independent
hash/receipt/entry-universe readback. Started full CA mapping after recording
resource estimates. Added a complete 526-taxon coverage table: 80 exact-sequence
reference proteins map to only nine Ascomycota taxa and 98 marker/taxon links,
including two hybrid taxa. Most coverage is concentrated in S. cerevisiae;
identical sequences across taxa do not provide independent experimental support.
This narrows the defensible generalization of the benchmark and motivates
broader homolog/reference evaluation. The prior turn completed direct geometry
comparisons and weighting summaries; this turn completes full acquisition and
expands mapping. Overall research scope remains active and incomplete.


### Experimental starting-model and chronology review completed

Added entry-level dependency flags for all 1,032 experimental references:
two explicitly list AlphaFold starting models, five other computational models,
548 only experimental starting models, and 477 have no starting-model annotation.
Preserved and verified all source lists and 2,843 entity/model date comparisons.
Target-chain attribution and training/template independence remain unresolved,
including for structures released after prediction creation. Full CA mapping is
confirmed running; phylogenetic stages still lack final completion receipts.
The preceding turn completed full coordinate acquisition and taxon coverage;
this turn contributes concrete evidence needed for circularity controls. The
complete project remains active and incomplete.


### Experimental-reference ESMFold control queue prepared and scheduled

Prepared all 45 experimental-reference proteins eligible for the current full-
length <=512-residue configuration; 35 longer proteins remain explicitly deferred.
Verified every emitted sequence and disjointness from four existing queues.
Estimated about 9.5 minutes of GPU time including the stated planning overhead,
then started a controller waiting for the pinned original GPU1 run to complete
cleanly and release the device. No extra GPU prediction process has launched.
The prior turn added experimental starting-model dependency flags; this turn
schedules the missing alternate-predictor evidence needed for three-way reference
comparisons. Full CA mapping and major phylogenetic/prediction jobs continue.
The overall goal remains active and incomplete.


### Experimental benchmark figure and coverage attrition verified

Generated and visually inspected a reproducible SVG/PNG/PDF figure showing
per-protein prediction–experiment agreement, deposition multiplicity and all
coverage exclusions. Every plotted protein median and count was recomputed
from comparison rows. The figure makes the change in eligible protein cohort
explicit (33 baseline, 32 at pLDDT 70, 18 at 90), preventing a confidence-filter
comparison from silently hiding reference attrition. Numerical measurements
remain descriptive pending experimental/context/training review. Existing
full CA mapping, species-tree and both GPU producers were verified live.
The previous turn scheduled missing experimental ESMFold controls; this turn
adds an inspected scientific figure and direct plot-data verification. The full
objective remains active and incomplete.


### Residue solvent-accessibility production started

Started a new analysis supporting the core/surface localization objective:
13,153 AlphaFold models and 6,910,765 residues, using 960-point Shrake–Rupley
with a 1.4 Å probe and four CPU workers. Prediction confidence and sequence
identity are retained. Isolated-chain accessibility is not treated as a validated
interface or biological core assignment. Analytic sphere/occlusion tests passed,
and the first eight completed tables passed sequence/grid/hash/area-sum readback.
The previous turn added the experimental benchmark figure; this turn starts
residue annotation across the full frozen AlphaFold marker snapshot. Full CA
mapping and other long-running stages continue. The full goal remains active.

### Local accessibility annotation and resolution sensitivity started

Extended the same residue-accessibility method to 5,121 frozen ESMFold models
and 1.18 million residues using two CPU workers. Eight initial local outputs
passed sequence/grid/hash/area-total checks. Started paired 960/3,840-point
numerical-resolution assessments on ten predetermined AlphaFold and six ESMFold
models, selected within fixed length strata independently of surface results.
The previous turn initiated full frozen AlphaFold accessibility; this turn adds
the second prediction source and numerical sensitivity. All stages remain
partial until their completion receipts and subsequent audits pass. The full
project goal remains active.


### Local solvent-accessibility resolution assessment completed

Completed and independently audited both sampling resolutions for six ESMFold
models and all 1,234 residues. Median absolute per-residue area difference is
0.320 Å², 95th percentile 1.137 Å², maximum 2.431 Å²; three residues switch
between exactly zero and nonzero area. These findings constrain categorical
burial interpretation without claiming biological accuracy or full-dataset
convergence. AlphaFold resolution sensitivity and both production annotation
runs continue. The previous turn launched local annotation and both sensitivity
runs; this turn finishes and validates the local sensitivity result. The full
research objective remains active and incomplete.


### AlphaFold accessibility-resolution assessment completed and audited

Waited on the verified live resolution process through the final long model,
then audited all ten selected models and 5,839 residues. Median absolute ASA
difference is 0.282 Å², 95th percentile 1.048 Å², maximum 2.049 Å²; 39 residues
switch between zero and nonzero area. Both source assessments are now complete,
with all outputs and predetermined model selections checked. This supports
continuous accessibility annotation while retaining explicit numerical and
biological uncertainty. Both full production annotation runs and the larger
research analyses remain incomplete. The prior turn completed the local
assessment; this turn completes its AlphaFold counterpart. The full goal remains
active.
