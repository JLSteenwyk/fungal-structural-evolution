# Progress and evidence

| Milestone | Current verified state | Evidence and remaining requirements |
|---|---|---|
| Repository and original objective | Public repository and reproducibility records maintained | Git history; docs/objective.txt |
| Sampling | 501 fungal entries plus 25 outgroups acquired; 21 incomplete fungal labels and two curated hybrids require identity-aware interpretation | metadata/analysis_manifest.tsv; docs/taxon-identity-sensitivities.md; species uniqueness not fully established |
| Assembly and protein QC | Broad QC for 526 and lineage-specific BUSCO for 501 fungi complete; FCS report inventory and exact coding-overlap audit complete | docs/assembly-quality-workflow.md; 518 usable reports, one checksum mismatch, seven external-source exceptions; biological review and omission/copy sensitivities pending |
| Species phylogeny and discordance | Both 526-taxon homogeneous guides audited; supported mixture analyses and gene trees running; 70/125 marker trees audited in latest frozen snapshot; supported-split guide conflict diagnostic computed | docs/phylogenetic-workflow.md; full support, discordance, root and sensitivity analyses pending |
| Families, domains and reconciliation | 64-taxon core complete; full 526-taxon artifact audit identifies one required family-tree repair, now running; full-proteome Pfam searches running | docs/orthology-workflow.md; docs/domain-annotation-workflow.md; reconciled families pending |
| Structural atlas | Frozen source cohorts include 13,153 AlphaFold and 10,048 combined ESMFold models; source-specific mapping, confidence and accessibility checks complete | docs/prediction-source-controls.md; larger acquisitions/predictions running; full atlas incomplete and source cohorts not pooled for inference |
| Sequence–structure analyses | Earlier 72-marker fits, conditional resampling and exploratory coupling complete; expanded AlphaFold 496 fits and benchmark audited; combined ESMFold 356 fits and corrected Gamma4 audit complete; 45,155 marker edges linked to both provisional guides with 180,620 audited branch estimates; FreeRate optimization/downstream comparisons running | docs/conditional-site-coupling.md; Earlier ESMFold 64 FCS omission fits audited; site-rate and coupling sensitivity pending; model, prediction and phylogenetic uncertainty remain |
| Coding-sequence analyses | 125 marker codon alignments, 1,655 nucleotide trees and 1,655 global MG94 diagnostics audited | docs/coding-sequence-workflow.md; all 13,240 nuisance-profile points audited; four optimization concerns and 29 FCS-exposed cases flagged; selection eligibility unresolved |
| Dating | Published summary chronograms inventoried and 27 calibration candidates catalogued | docs/dating-workflow.md; specimen/placement/prior review and joint age uncertainty pending; no time-normalized project rates |
| Evolutionary objectives 1–8 | Intermediate analyses available; no objective set declared complete | docs/research-plan.md; final branch/clade tests, duplication/domain/ecological analyses, selection and ancestral case studies remain |
| Figures, methods and mechanistic cases | QC and exploratory figures and methods available | docs/figures/; docs/methods-draft.md; final integrated results and validation proposals incomplete |

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


### Second TFIIB repeat supports the fixed class-associated split

Repeat 2 completed with all 1,040 gene-copy tips and 1,000 retained bootstrap
trees. The fixed combined-domain bipartition is present with SH-aLRT 99.8 and
998/1,000 exact bootstrap occurrences. An independently checked named-tip search
confirms the ML split. A nearby split improving annotation agreement by one
Vairimorpha necatrix copy is much weaker (65.7 SH-aLRT, 471/1,000 bootstraps),
so both results are retained explicitly. The result reinforces the copy-mixture
review without certifying orthology or independent replication. Repeat 1 and
species/gene reconciliation remain pending. The prior turn completed AlphaFold
surface-resolution validation; this turn advances phylogenetic sensitivity.
The full project remains active and incomplete.

### Frozen structural similarity clustering started

Started Foldseek clustering across all 18,274 frozen AlphaFold/ESMFold models,
retaining distinct source/model identities and original taxon/marker mappings.
Verified all aliases, symlink targets and source hashes. Explicit alignment,
coverage, TM-score, confidence-seeding and clustering settings are pinned with
resource estimates; database construction is running. Similarity groups remain
exploratory and require representative/member and threshold validation, and are
never equated with orthogroups. Updated the top-level atlas milestone to reflect
current evidence rather than the early 128-model checkpoint. The preceding turn
validated the second TFIIB repeat; this turn initiates the frozen-model atlas
clustering stage. Full CA mapping and other long jobs continue; the complete
research goal remains active.


### Both TFIIB repeats and initial structural clustering completed

Finished the repeat-1 fixed-split assessment: 998/1,000 exact bootstrap
occurrences, matching repeat 2, with the same 485/555-tip ML bipartition confirmed
by an independent named-tip search. The combined alignment has 1,000/1,000;
repeat results are sensitivity checks, not independent replication. Reconciliation
and the mixed-copy marker restriction remain required.

The 18,274-model clustering also completed: 2,249 groups, 1,276 singletons,
302 mixed-source groups. All database identities and final membership/provenance
fields passed complete readback. Representative/member alignment criteria and
biological interpretation remain unvalidated. The prior turn launched clustering;
this turn completes its first result and the remaining repeat sensitivity.
The full research goal remains active and incomplete.

### Structural cluster edge validation finds memberships requiring review

Realigned all 32,050 directed representative/member edges; all returned.
15,043 pairs pass the reported criteria bidirectionally, 979 fail both directions
and three pass one direction. The discrepant memberships affect 200 original
groups. Readback also found 102 approximate alignment-TM scores above one,
explicitly flagged rather than clamped. Started exact-score conversion of the
retained alignments before conservative membership decisions. The prior turn
completed clustering and membership identity checks; this turn tests actual
edge criteria and exposes limits that matter for atlas interpretation. Original
clusters remain exploratory, and the full research goal remains active.

### Structural groups linked to complete frozen taxon/marker provenance

Attached and independently verified all 18,815 model–taxon–marker links across
2,249 original similarity groups. Twenty-three groups have multiple marker labels,
95 include outgroups, and all 200 groups with approximate edge failures retain
explicit flags. This creates a traceable connection to the phylogenetic inputs
without treating clusters as orthogroups. Exact-score conversion is verified
live; its downstream reviewer is prepared but awaits completion. The preceding
turn tested edge criteria; this turn completes group provenance annotation while
preserving unresolved validation. The full research goal remains active.

### Full experimental residue mapping audited; full agreement run started

Completed all 1,032 entry mappings and audited 1,306,701 exported position
records, including 1,008,964 unambiguous full-occupancy CA positions. Five
selected raw mmCIF entries passed independent field readback; the sole
zero-occupancy exclusion was separately confirmed in 5LQW. Started full
prediction–experiment geometry after the audit passed, retaining the same
protocol and independent numerical recalculations as the partial run.
Experimental quality/context/training independence and broad taxonomic coverage
remain unresolved. The previous turn annotated structural groups; this turn
completes the full reference correspondence dataset and expands its benchmark.
The overall project goal remains active and incomplete.


### Full experimental agreement and exact cluster-score review completed

Completed and independently recalculated all 8,950 accepted experimental
comparisons. Published protein-weighted and fixed-cohort summaries plus the
full reference figure (72 baseline-eligible proteins of 80). Deposition
multiplicity materially changes descriptive medians; limited reference coverage
and training independence remain unresolved. Generalized plotting labels and
corrected inherited partial-snapshot summary wording without altering old runs.

Exact cluster-score conversion also completed: 14,769 representative/member
pairs pass both directions; 1,256 do not. Seventy-two directed alignment-TM
scores remain above one and explicitly fail review, with their values preserved
for investigation. Original groups remain exploratory. Existing ESMFold,
phylogeny and orthology jobs were verified live. The previous turn launched
full agreement; this turn completes its summaries and exact edge review.
The overall research goal remains active and incomplete.

### Reviewed representative groups derived and audited

Derived conservative memberships from the completed exact edge review: 17,018
models retained, 1,256 deferred, with all 18,815 taxon–marker links preserved.
Every retained nonself member passes both directed comparisons to its original
representative. Deferred models receive no derived family assignment. Complete
readback recalculated eligibility from raw scores and verified model provenance,
partitions and annotation links. This advances atlas usability without equating
similarity groups with orthogroups. The prior turn completed exact score review;
this turn makes its eligibility decisions explicit in downstream group tables.
Existing prediction, guide-tree and orthology processes remain verified live.
The full research goal remains active and incomplete.

### Structural filtering coverage quantified across the full manifest

Summarized all 526 taxa, 125 marker labels, available lineages and both prediction
sources before and after representative-edge filtering. All 459 represented
taxa retain models; 67 have no frozen input model. Retention differs by source
(94.6% AlphaFold, 89.4% ESMFold), with representatives reported separately from
tested members and source/length/confidence confounding explicit. Independent
readback verified all 677 taxon, marker, lineage and source count rows. The
preceding turn derived reviewed groups; this turn characterizes their missingness
for downstream evolutionary interpretation. Prediction, paired-fit, accessibility
and paired-geometry processes were observed live; their full receipts are not
yet available. The overall research goal remains active and incomplete.

### Pinned source and full alignment traces support normalization explanation

Reviewed Foldseek source at the reported binary revision and exported all 32,050
saved alignment traces. Every CIGAR agrees with inclusive endpoints, whereas
alignment-TM normalization uses end-minus-start. All 72 scores above one have
fewer denominator positions than matched positions and obey the resulting
upper bound. This provides a concrete explanation to test with corrected
rescoring; no scores were clamped, rescaled or reinstated. The prior turn
quantified filtering coverage; this turn identifies a likely implementation
cause behind some exclusions. Existing prediction, paired-fit and geometry
processes remain live. The overall research goal remains active and incomplete.

### Matched-toolchain normalization control build launched

Started an isolated four-job build of the pinned Foldseek revision, producing
an unmodified control and an inclusive-span output-normalization variant with
the same toolchain. Source/configuration hashes and a resource plan are tracked;
configuration passed and the unmodified compilation is verified live. Prepared
receipt-gated rescoring of the unchanged saved alignments. The previous turn
identified a likely normalization issue; this turn starts the controlled test
needed before correcting scores or membership. Existing prediction and paired
analysis processes remain live. The overall research goal remains incomplete.

### Normalization-control builds completed and full rescoring started

Both pinned-source builds completed. Binary/source/patch readback passed, and
the production executable remains unchanged. Launched two sequential exact
conversions on all saved alignments, with the unmodified control verified live.
Prepared a reviewer that separates toolchain effects from normalization effects
and rejects changes to other score fields between matched builds. The preceding
turn launched compilation; this turn completes the build and starts the actual
controlled score comparison. Results and corrected memberships remain pending.
The full research goal remains active and incomplete.

### Full frozen local-prediction domain comparisons launched

Started conserved-domain geometry across 266,788 baseline pairs from all 5,121
frozen ESMFold models, using completed Pfam annotations and local PAE exports.
The existing three geometry tests pass, and source receipts/resource estimates
are pinned. This extends the earlier small AlphaFold domain analysis to the full
available local-source collection, separating within-domain geometry from
whole-protein orientation effects. The producer and normalization-control
conversion are both verified live. The previous turn started controlled score
rescoring; this turn advances an independent prerequisite for interpreting
structural change. Full domain outputs and the overall research goal remain
incomplete.

### Expanded structural-feature dependency audits started; local audit complete

Completed the feature-overlap audit of all 72 frozen ESMFold paired markers:
4,669 taxon–marker alignments and 714,936 observed features. Of 3,901,410 pairs
sharing coordinates, 44.24% exceed the reach of a length-10 circular block and
25.75% that of length 30. All paired alignment identities and observed counts
passed FASTA readback. Expanded AlphaFold feature auditing, domain geometry and
normalization controls remain verified live. The preceding turn extended domain
geometry; this turn advances uncertainty assessment for sequence–structure
comparisons. Potential dependence is not measured covariance, and the full
research goal remains active and incomplete.

### Expanded AlphaFold dependencies completed and both source counts verified

Completed the 124-marker AlphaFold feature audit: 13,565 taxon–marker alignments,
4,148,852 observed states and 23,393,290 pairs sharing coordinates. Fractions
beyond length-10 and length-30 circular blocks are 46.81% and 28.97%. Added a
reusable complete FASTA/count/fraction readback and ran it successfully on both
expanded sources. This completes the dependency measurement started last turn;
it does not establish calibrated branch intervals or statistical covariance.
Domain comparisons and normalization controls remain verified live. The full
research goal remains active and incomplete.

### Unmodified normalization control completed with identical reported output

After a verified wait, the unmodified control completed all 32,050 directed
alignments. Every exported numeric field matches the production exact-score
conversion, ruling out observed toolchain differences at reported precision.
The patched conversion has started and remains live; domain geometry also
continues. This advances beyond the preceding verified-wait turn by completing
and checking the first controlled result. The overall goal remains active and
incomplete; no corrected score or membership claim has been made.

### Full frozen ESMFold domain run completed with full candidate accounting

Completed 171,474 conserved-domain comparisons over 135,959 marker–taxon pairs
and 58 Pfam domains. All 266,788 baseline pairs and 306,789 candidate dispositions
are accounted for: 128,146 no-shared-eligible-domain exclusions and 7,169
insufficient-coverage domain exclusions. Added and executed full correspondence,
count, baseline-field and fit-relation readback on both the full local dataset
and earlier comparison dataset. Independent coordinate checks and interdomain
placement review remain pending. The previous turn completed the unmodified
Foldseek control; this turn completes a substantive domain-level dataset.
The overall research goal remains active and incomplete.

### Independent domain geometry verified; full interdomain analysis started

All 58 identity-selected domain comparisons passed independent SciPy rotation
and distance recalculation, one per Pfam accession. Started placement/confidence
analysis of all 35,708 eligible domain-pair/taxon combinations at three PAE
thresholds. This follows the previous turn's completed domain dataset and
complete candidate accounting with a separate geometry implementation and
an analysis of interdomain uncertainty. Sampled checks do not validate every
row, and interdomain results remain pending. The overall research goal remains
active and incomplete.

### Full interdomain confidence completed; local-domain figure published

Completed all 107,124 interdomain rows and verified every expected domain-pair
combination at all three PAE thresholds. At PAE10, the median retained
cross-domain fraction is 0.825 and 195 combinations retain no confident pairs.
Published a visually inspected figure of all 171,474 domain comparisons and
its descriptive summaries. These distinguish within-domain fit from placement
uncertainty without claiming evolutionary rates or biological flexibility.
The previous turn launched placement analysis; this turn completes its result
and full grid readback. The overall research goal remains active and incomplete.

### Patched normalization completes and resolves out-of-range scores

After verified waits, the patched conversion completed and passed full
32,050-alignment control review. The unmodified build matches production;
only alignment-TM changes under the patch. All 72 out-of-range values disappear,
while 68 directed rows change fail→pass and 46 pass→fail. Membership revisions
still require bidirectional review. The previous turn verified a live wait;
this turn completes the controlled normalization result. The full research
goal remains active and incomplete.

### Corrected scores propagated into fully checked structural memberships

Reviewed both directions of all original edges and regenerated groups:
17,029 retained models and 1,245 deferred. Complete raw-score and provenance
readback passed; 34 restorations and 23 new exclusions match direct membership
set differences. Updated coverage retains all 459 represented taxa, with absent
inputs explicit in the full manifest. Earlier versions remain immutable.
The previous turn verified the normalization correction; this turn propagates
it into usable group and coverage tables without treating groups as orthology.
The overall research goal remains active and incomplete.

### All local paired fits completed; full resampling launched

Completed and audited 288 fits across 72 local-source markers and 9,122 paired
branches. Warnings and near-zero branch counts remain explicit (248 and 208
fits, respectively). Started 43,200 paired resampling draws over every marker,
with eight workers, source/configuration pins and a recorded resource estimate;
all three resampling tests passed. The prior turn updated corrected structural
groups; this turn advances branch estimation and its conditional uncertainty
analysis. Nonlocal feature dependence and other uncertainty limits remain,
and the overall research goal stays active and incomplete.

### Complete local tree-path/geometry point benchmark assembled

Linked all 259,780 accepted and 6,813 excluded geometry pairs to sequence and
three structural-model tree paths across 72 markers. Full pair accounting,
tree provenance/topology checks, 2,848 independent path traversals and complete
original-field preservation readback passed. This follows the preceding turn's
completed local fits and active resampling with the direct structural benchmark
needed for later coupling analysis. Point estimates remain separate from
uncertainty, additivity and phylogenetically adjusted inference. The overall
research goal remains active and incomplete.


### Local tree-path benchmark summarized across all 72 markers

Computed 864 descriptive within-marker rank associations for four path models
and three geometry metrics. All correlations and cohort counts passed separate
SciPy recalculation, and all 12 equal-marker summaries passed readback. Local
geometry has stronger median rank association than whole-protein RMSD across
these models; ancestry, uncertainty and confidence effects remain unadjusted.
Both GPU prediction jobs, expanded AF fits/geometry and full local resampling
were verified live. This advances the completed point benchmark; phylogenetic
coupling tests and the overall research goal remain incomplete and active.


### Accessibility annotations checked against original atom records

Added a reusable full-production audit with an explicit frozen-partial mode.
All 4,652 entries available at launch passed raw mmCIF sequence, residue grid,
atom-count, CA confidence, hash and ASA-total checks: 1,059,944 residues. This
prepares annotations for subsequent structural localization analyses; it does
not recalculate ASA or validate biological core/surface assignments. The 469
entries outside the frozen audit still require checking. Larger phylogenetic
and prediction runs were verified live. The previous turn completed the local
rank benchmark; this turn advances residue annotation integrity. The overall
research goal remains active and incomplete.


### Full local benchmark visualized without marker subsampling

Produced a three-panel figure showing every one of 864 marker/model/geometry
correlations and the 12 equal-marker medians. Paired lines preserve marker
identity across four models. All median values and SVG point counts passed
readback, and the PNG was visually inspected. This follows the prior residue
annotation audit and advances the reproducible benchmark figure deliverable.
Larger geometry, fit, prediction and resampling processes remain live; the
research goal remains active and incomplete.


### Full accessibility audit queued behind verified live production

After the preceding verified wait, production was confirmed live at PID 607025
with 4,895 completed local-source entries. Added and launched a controller that
waits for this specific producer, requires its final receipt and pinned inputs,
and then audits all 5,121 models without partial mode. Runtime and resource
estimates are recorded; the existing host is used. The controller is live and
waiting, so full-audit completion is not yet claimed. This automates the next
necessary validation step while broader analyses continue. The overall research
goal remains active and incomplete.


### Full local accessibility completed and linked to all paired sites

Following the preceding verified wait, all 5,121 models and 1,183,341 residues
completed production and passed the queued full raw-coordinate audit. A new
projection links all 714,936 observed paired AA/3Di sites across 72 markers,
4,669 taxon/marker cells and 4,630 models. Full alignment-site set and annotation
readback passed, and partial-audit inputs were confirmed to fail before output.
This advances the evolutionary localization dataset; normalized exposure,
core/surface sensitivity and phylogenetic tests remain unfinished. The overall
research goal remains active and incomplete.


### Exposure normalization added with explicit scale sensitivity

The preceding turn completed full local accessibility and paired-site projection.
This turn normalized every one of 714,936 sites with two published reference
scales, retaining original annotations and unclipped indices. Full arithmetic,
source-field and summary readback passed. Between 1.16% and 8.90% of sites change
sides of the five diagnostic thresholds across scales; no biological core/surface
labels or evolutionary transitions are inferred. DSSP/ShrakeRupley equivalence
remains unestablished. Main fit, geometry, GPU prediction and resampling jobs
were verified live. The overall research goal remains active and incomplete.


### Every paired site linked to minimum tree-based character changes

Following complete exposure normalization, added a topology-aware diagnostic
for all 16,571 columns and both alphabets across 72 local-source markers.
Minimum-change scores and extant exposure summaries preserve all 714,936
observations. Exhaustive small-tree tests passed, and every one of 33,142 real
scores matched a separate recurrence; all site exposure summaries also matched.
These are minimum counts with explicit model limits, not branch rates or
phylogenetically adjusted coupling tests. Main jobs were verified live. The
overall research goal remains active and incomplete.


### Live resampling missing-data event reconstructed

After the preceding completed site-parsimony dataset, verified all main jobs
live and inspected a newly reported unestimable local resampling draw. Complete
AA/3Di reconstruction confirms that block sampling leaves F456999 entirely
missing despite 94 observed original sites. The intended no-fit/no-replacement
handling is intact. Resampling remains active and requires a final all-draw
summary. The overall research goal remains active and incomplete.


### Site topology sensitivity launched over every saved AA bootstrap tree

Following the preceding verified wait, launched original-character minimum-change
calculations over all 72,000 local AA bootstrap topologies, retaining all score
arrays. Eight workers and resource forecasts are recorded. The process was
verified live with 14 marker outputs already complete; full-stage completion
and independent output auditing remain pending. This adds topology sensitivity
to the exposure-linked site dataset alongside ongoing column resampling. The
overall research goal remains active and incomplete.


### Full topology sensitivity completed and validated

Following the preceding launch, all 72,000 AA bootstrap topologies completed.
The 33,142,000 stored scores passed full grid/bound/provenance and summary checks;
one predetermined tree per marker supplied 33,142 independently recalculated
scores, all matching. Site minimum counts vary across topology samples at 7,514
AA and 4,965 3Di sites. These remain descriptive topology sensitivity, not
calibrated confidence or acceleration claims. Column resampling and other major
stages remain live. The overall research goal remains active and incomplete.


### Conditional site-rate estimates completed across all local markers

After the preceding verified wait, completed 288 fixed-topology empirical-Bayes
rate fits and exported 66,284 site estimates. Corrected a parser mismatch with
actual IQ-TREE3 column labels in a new immutable run; both schema and rejection
tests passed. All source/model/topology identities, rate grids/category bounds
and likelihood accounting passed audit. Warnings in 248 fits remain explicit.
This adds a model-based evolutionary measure for later exposure/coupling analyses;
model adequacy and inferential controls remain pending. The overall goal stays
active and incomplete.


### Matched site-evolution table assembled and fully checked

After the preceding verified wait, joined all 16,571 sites to 66,284 rate
estimates, exposure/parsimony/topology diagnostics and confidence/coverage/model
covariates. Every inherited and rate field passed readback; all new confidence,
model and marker summaries were independently checked. This prepares consistent
inputs for controlled coupling/exposure models without treating the assembly
as statistical adjustment. Main long-running stages remain live. The overall
research goal remains active and incomplete.


### Symmetry screening completed with unavailable results explicit

Following the preceding matched analysis table, screened all 144 paired local
alignments without removing markers. Independent reconstruction of all 533,186
Bowker pair counts matched the tool outputs. Of 432 maximum-test results, 186
are unavailable; only 7/72 3Di alignments have estimable marginal/internal
maximum tests. No model-adequacy conclusion follows from unavailable or nominal
results. Calibration, multiplicity and further model sensitivities remain
pending. Main long-running stages remain live; the overall goal stays active
and incomplete.


### Full local resampling audit queued automatically

The 72-marker, 43,200-draw producer remains live. A controller now waits for
its exact process identity to exit, requires its complete execution receipt,
checks pinned configuration and code, and runs the full existing resampling
auditor into a new immutable output. It does not restart or modify the producer.
The audit reconstructs both alphabets from the sampled columns, checks artifact
hashes and fixed-topology branch estimates, and retains unestimable draws.
One CPU worker, 8 GB memory and 1 GB output are planned, with a conservative
1–12 hour audit estimate on the existing host. Configuration is recorded in
`metadata/esmfold_paired_resampling_audit_controller_config.json`.
The controller compiled and entered the expected waiting state. The audit has
not run yet; interval calibration and the overall research goal remain incomplete.


### Full local rate-heterogeneity sensitivity launched

Started all 288 FreeRate4 fits on the original paired alignments and fixed AA
topologies. This extends the model sensitivity needed before interpreting
sequence–structure coupling. Initial fits succeeded; full execution, FreeRate
audit and matched comparison remain pending. The revised auditor passed all
288 existing Gamma fits with identical prior statistics and site/warning tables.
Prepared the matched site/branch comparison with rank and model-identity checks.
Resources, provenance and reproducible commands are documented in
`docs/site-specific-evolutionary-rates.md`. Main long-running stages and the
resampling audit controller remain live. The research goal stays incomplete.


### Functional correspondence added to the matched site-evolution frame

Joined existing annotations to all 16,571 local sites without changing prior
estimates. The 1,189 observed annotation links represent 1,066 taxon-site cells
at 28 sites; 10 sites across seven markers carry conserved candidates. Full
field and count readback passed. Unannotated sites remain unknown, and repeated
profile links are not counted as separate taxon observations or events. This
prepares focused functional-site case studies while exposing limited annotation
coverage. FreeRate and other major jobs remain running; the goal is incomplete.


### FreeRate comparison completed; seven optimization diagnostics retained

After verified waits, all 288 FreeRate fits completed. The full audit passed
66,284 rate rows and retained warnings from 249 fits. All 66,284 matched site
and 36,488 matched branch values passed readback, as did all rank correlations,
median differences, tree totals and likelihood changes. Median rank agreement
is high across model specifications, but seven FreeRate fits have lower
likelihoods than Gamma (largest decrease 17.3648 log units). These remain
explicit unresolved optimization diagnostics; no likelihood-based selection,
acceleration or model-adequacy conclusion is made. Other major jobs and the
full research goal remain incomplete.


### Lower FreeRate likelihoods traced to optimization sensitivity

Reviewed the pinned IQ-TREE3 source and ran 28 explicit-start diagnostic fits
for all seven flagged cases. All outputs and 6,348 rate entries passed full
audit. Every case has a diagnostic likelihood above Gamma; the original worst
case improved by 38.1536 log units. This supports initialization/optimizer
sensitivity, not a biological reversal. Original outputs remain preserved and
the diagnostic best-fit table is separate. Consistent full-run optimization,
updated sensitivity inference and the broader research goal remain outstanding.


### Consistent optimization launched across all local FreeRate fits

Extended the exact four-refit procedure to every one of the 288 model fits,
creating 1,152 requests. The full request grid and retained settings passed
readback; all 28 earlier diagnostic requests match except output location.
Initial fits are completing. Full execution, output audit and revised model
sensitivity remain pending. Resampling is also still live with 209/216 batches
complete at the launch check. The overall goal remains active and incomplete.


### All full-local paired resampling draws complete; full audit started

The 72-marker run completed all 216 marker/block batches and 43,200 planned
paired draws: 43,184 reported estimable draws and 16 retained unestimable draws.
The full marker/block grid, every batch receipt hash, replicate receipt-name
grid and accounting totals passed readback. The producer exited successfully.
Its controller then launched the full raw-output auditor automatically; live
process and initial audited batches were verified. Full audit and conditional
interval summaries remain pending. These completion counts alone do not
validate branch uncertainty or establish biological acceleration. Other major
production stages and the broader goal remain incomplete.


### Joint tree-path uncertainty prepared behind the full audit

While the raw resampling audit progresses, implemented propagation of joint
branch draws into all accepted/excluded exact-site tree paths at block lengths
1/10/30. Three covariance/validation tests passed; the missing-audit gate
correctly prevented full-data execution and created no output. Resource plan
and reproducible commands are recorded. Full audit, path execution and readback
remain pending, as do the broader evolutionary analyses.


### Automatic handoff queued from resampling audit to joint path summaries

The full resampling audit remains live. A pinned controller now waits for that
specific auditor process, requires its complete receipt, and launches the tested
joint-path summary script automatically. The controller compiled and entered
the expected waiting state; no dependent output was claimed complete. Full path
execution/readback and the wider research goal remain outstanding.


### Full resampling audit and joint path readback completed

The full local resampling audit passed all 86,368 fitted trees from 43,200
paired draws (16 unestimable) and wrote 54,732 branch interval rows. Warnings
remain recorded for 74,465 fits. The queued successor completed joint path
summaries for all 259,780 accepted pairs and 6,813 geometry exclusions across
72 markers and block lengths 1/10/30. All source batch pins and the configuration
were rechecked against the resampling receipt.

A new reproducible readback verified every inherited pair field, cohort identity,
summary count, interval bound and covariance bound. For one SHA256-selected pair
per marker, independent Bio.Phylo traversals of every estimable draw reproduced
the percentiles, standard deviations and paired covariance: 216 pair/block
checks and 86,368 tree traversals passed. Numerical summaries for all other
pairs were not independently recomputed. The first readback invocation hit a
Python-version incompatibility in hashlib.file_digest before any readback output
was created; streaming SHA256 fixed compatibility and the full rerun passed.

Receipts and methods are tracked with the script. These remain conditional
sampling sensitivities, not calibrated confidence intervals or evolutionary
acceleration results. The consistent 1,152-fit FreeRate optimization remains
live, alongside the broader phylogeny, orthology, domain annotation, prediction
and structural analysis jobs. Full project completion remains outstanding.


### Full optimized model-comparison handoff queued

All major production PIDs were revalidated live. Prepared the matched-rate
comparison to select the best observed R4 fit consistently across all 288
marker/model combinations, preserving original fits when better or tied.
Default-mode full-data compatibility passed with all three tables byte-identical;
the seven-case partial diagnostic input was explicitly rejected without output.
Compiled and launched a pinned controller that waits for the exact live full
optimization process, then audits all 1,152 refits and runs the revised
comparison. The controller entered its expected waiting state. Full execution,
selected-fit numerical readback and broader evolutionary analyses remain
outstanding; the overall goal is active and incomplete.


### Full raw rate-comparison readback implemented and validated

Added a reproducible checker for all source rates, tree edges, likelihoods,
summary statistics and fit choices. The original complete comparison passed
all 288 fit summaries, 66,284 site comparisons and 36,488 branch comparisons.
The checker uses separate rate parsing and postorder edge collection, with
shared Bio.Phylo Newick parsing explicitly recorded. Optimized fit choices
will be checked from all four raw diagnostic reports rather than trusting the
selected-fit summary. Full optimization and its queued audit/comparison remain
live; execution of the optimized readback remains pending. No biological
rate conclusion or project completion is claimed.


### Full 513–768-residue marker queue prepared and checked

Prepared all 5,512 canonical sequences in this length band from the original
and follow-on candidate queues, preserving 5,549 taxon/marker links across
279 taxa and 108 markers. Updated exact-sequence reuse found two available
models, leaving 5,510 candidates. Full sequence/link/reuse and resource-scenario
readback passed. Frozen timing observations yield explicitly extrapolated
41.82–66.70 GPU-hour inference scenarios and a 100.05-hour scheduling allowance.
No longer-protein prediction was launched or queued yet. Existing GPU jobs,
full FreeRate optimization and its queued audit remain live. Larger deferred
proteins, broader structural coverage and the full research goal remain
incomplete.


### Longer-marker GPU handoff queued

Compiled and launched a pinned controller for the full 513–768-residue queue,
waiting behind the already queued experimental controls on GPU 1. The exact
predecessor process is live and the new controller entered its expected waiting
state. It requires clean predecessor completion, refreshes model reuse and
resource observations into a new immutable input snapshot, verifies the complete
selected sequence/link universe, then launches the remaining candidates with
max-length 768 and the existing OOM stop. The initial queue has 5,510 pending
proteins; launch-time reuse can change that count. No longer-protein inference
has started yet. Existing prediction and model-fitting jobs continue; the
full project remains incomplete.


Full lineage BUSCO execution and raw-output audit completed: 501 fungal entries,
998,331 marker calls and 943,724 protein-hit rows. The versioned quality table
retains all 526 entries and outgroup broad-panel scores. All exported artifact
checksums match the audit receipt. Six panel-specific summaries preserve their
distinct denominators; no universal cutoff or biological-loss claim is applied.
See docs/lineage-completeness-workflow.md and metadata/lineage_busco_audit_receipt.json.


Full FreeRate optimization, automatic audit/comparison and independent raw-output
readback completed: 1,152 refits across all 288 source fits; final comparisons
contain 66,284 sites and 36,488 branches. All 288 selected refits improve the
original FreeRate likelihood and remove the earlier seven >0.1 Gamma deficits.
Model sensitivity remains conditional; this is not biological model adequacy or
an acceleration result. See docs/site-specific-evolutionary-rates.md and
metadata/esmfold_optimized_rate_comparison_readback.json.


Queued full original/follow-on prediction artifact audits and mmCIF conversion,
with exact live producer identity and code/input pins. Expected full cohorts are
10,522 and 4,252 models. The strengthened auditor's seven tests passed, including
metadata corruption and completion accounting checks. Both handoff controllers
were observed live waiting for their still-running producers; all 11 dependency
pins per configuration matched. Full audit/conversion results remain pending.
See docs/structure-prediction-workflow.md and the two snapshot controller configs.


Prepared and verified full-matrix PMSF resource scenarios from both current guide
startup reports. Both 526-tip matrix manifests and taxon grids passed readback.
C20 planning memory plus 25% allowance is approximately 417k/542k reported MB;
C60 can exceed host RAM. A v1 512G allowance was corrected to 600G in v2 before
launch. Proposed full-matrix crossed-guide sensitivity retains all taxa and
requires completed guide receipts, serial execution and 750 GiB available memory.
Resource evidence and limitations are in metadata/species_mixture_resource_receipt.json
and docs/phylogenetic-workflow.md. This stage remains unlaunched; existing guide
and gene-tree producers were confirmed live.


Follow-on ESMFold inference completed: 4,252/4,252 eligible sequences, no OOM or
interruption; producer session 43888 exited zero. Full artifact readback started
(PID 3564443 under snapshot controller 3228696/session 77902), and the ecology
controller launched the authorized 675-sequence queue (inference PID 3566770).
Both launch identities were observed live and the ecology predecessor hash
matched the completed chunk. Audit/conversion and ecology results remain pending.
Versioned execution/launch receipts are linked in docs/structure-prediction-workflow.md.


Profile guide completed and passed full 526-tip/report-edge audit; 520 nominal
composition failures retained. First full-matrix C20-PMSF sensitivity launched
with 1,000 SH-aLRT and 1,000 UFB replicates, 16 threads and 600G limit
(parent PID 3700451, IQ-TREE PID 3700949, session 91598). Initialization estimates
333,579 MB. Other crossed-guide runs and full supported-tree audit remain pending.
Follow-on model audit also completed: 4,252 models, 4,304 links, 84 taxa, 98 markers;
all exported identities/links checked independently. Conversion is live under
PID 3676262/controller 3228696. See the phylogenetic and prediction workflows.


Queued the full 675-model ecology audit/conversion with explicit ecology marker
links (controller session 55711). All candidate identities, input hashes and link
coverage were checked; the live producer is pinned by PID/start time. The added
controller differs from the existing pinned controller only by explicit link
forwarding and its description. Original/follow-on live controller files remain
unchanged. Follow-on conversion reached 1,500 models; the full-matrix PMSF job
was confirmed live. Completion of these stages remains pending.


Queued full follow-on residue mapping after successful conversion, preserving
all audited model and originating marker-link identities. The controller binds
all 125 profile receipts to the completed alignment collection; all 134 source
pins passed readback. Controller session 78949 is live waiting for conversion
controller 3228696. Mapping, residue readback and PAE/native qualification remain
pending. See metadata/esmfold_followon_mapping_controller_config.json and the
structure-prediction workflow.


Full follow-on conversion completed and controller session 77902 exited zero.
All 4,252 inventory/provenance identities, converted-coordinate hashes and original
prediction-receipt hashes passed independent readback. The residue mapper launched
successfully (PID 175824 under controller 4137875/session 78949). Mapping and
subsequent residue/PAE/native-feature checks remain pending. Completion and launch
receipts are versioned in metadata/esmfold_followon_*.


Full follow-on residue mapping and independent readback completed: 4,252 models,
4,304 marker links, 84 taxa and 1,017,794 residues. The new reusable auditor also
passed the entire prior local snapshot (5,121 models, 900,303 residues). All
expected retained positions, matrix amino acids, original NPZ confidence and
link/model identities checked. Mapping completed during queue setup, so the
prelaunch wait attempt exited and direct readback proceeded after terminal
verification; no producer was restarted. PAE/native qualification remains pending.
Receipts are metadata/esmfold_followon_full_residue_readback.json and
metadata/esmfold_existing_full_residue_readback.json.


Full follow-on native extraction completed (session 24784 exit zero); executable
and all four native source hashes match the prior validated encoder. Directional
PAE export is live (session 31256), and coordinate-feature reconstruction/audit
started (session 91321), covering all 4,252 models and 1,342,046 residues. PAE
export covers 476,502,626 directional entries. Joint PAE/native qualification
and updated paired inputs remain pending. Resource plans and commands are recorded
in the prediction workflow and metadata/esmfold_followon_*.


Expanded AlphaFold direct paired-site geometry completed for all 124 ready markers:
737,851 accepted and 7,002 excluded pairs (744,853 total). Producer session 73482
exited zero. Full pair-grid/character audit launched under PID 812449/session 24060,
with one SHA-selected numerical geometry/PAE check per accepted marker. Scope and
resource limits are recorded in docs/tree-path-geometry.md; full numerical
recomputation of every geometry is not claimed. Geometry audit and fitted-tree
comparisons remain pending. Follow-on coordinate/PAE validation continues live.


Expanded AlphaFold paired-site geometry audit completed (session 24060 exit zero):
all 744,853 pair records checked; 737,851 accepted and 7,002 excluded. Independent
SciPy geometry and directional PAE indexing matched for 124 preselected pairs,
one per accepted marker. All source receipt hashes and the 124-marker sample grid
were rechecked. Scope remains full pair/character verification with sampled
numerical geometry, not full numerical reoptimization or branch-rate inference.
Receipt: metadata/gdm_expanded_paired_geometry_audit_receipt.json.


Full follow-on coordinate feature audit completed (session 91321 exit zero):
4,252 models, 1,342,046 residues, 1,333,542 valid states and 8,504 invalid terminal
states. Independent exported-encoding readback checked all 1,333,542 six-residue
confidence minima and every model/summary count. Joint PAE qualification is queued
under session 72300 after exact exporter PID 642606; all ten dependency pins
matched. Pre-PAE six-residue pLDDT coverage is 1,008,659 states, not final qualified
coverage. PAE export, qualification readback and paired inputs remain pending.


Full expanded AlphaFold accessibility calculation completed: 13,153 models and
6,910,765 residues; producer session 64463 exited zero. Compact execution receipt
preserves the full 13,153-entry receipt-map checksum. Started the full raw-mmCIF
and residue-table audit (session 22787); independent ASA integration is not part
of this readback. Exact paired-site projection/normalization remains pending.
See docs/residue-accessibility.md and metadata/gdm_full_accessibility_*.

Queued the expanded AlphaFold paired-site accessibility projection and reference
normalization behind exact audit PID 1312606 (controller PID 1888349, session
20643). The controller pins 16 dependencies, requires full matching audit counts,
checks resource headroom and preserves immutable outputs. Syntax/import checks
and the live wait state passed; neither downstream calculation has run yet.
Audit progress reached 3,500/13,153 models; follow-on PAE export reached
3,000/4,252 models. This handoff prepares exposure covariates for all 124 ready
AlphaFold markers; controlled evolutionary analyses and independent output
readback remain pending.


Started full follow-on ESMFold accessibility for 4,252 models and 1,342,046
residues (PID 2004958, session 50918), with the same 960-point method as earlier
cohorts. Resource estimates use measured earlier-cohort worker times, with wider
planning allowances. The established full-output auditor is queued under
session 61482 after exact producer exit; 11 dependency pins are recorded.
Initial models produced complete residue tables. Neither full calculation nor
full audit is complete. See docs/residue-accessibility.md and
metadata/esmfold_followon_accessibility_*.


Follow-on PAE export completed for all 4,252 models, zero failures (session 31256
exit zero). Full manifest identity/source linkage and all 3.81 GB compressed-file
hashes passed a separate readback. Exact matrix round-trip was checked by the
exporter, not repeated independently. Qualification automatically launched as
PID 2283888 under controller session 72300 and remains live. Final confidence
counts, independent context readback and paired-input preparation remain pending.
See docs/structure-prediction-workflow.md and metadata/esmfold_followon_*pae*.


Follow-on qualification finished: all 4,252 models, 1,006,992 full-chain residues
passing six-residue pLDDT70 and PAE10. Independent original-NPZ context readback
passed all 1,333,542 valid follow-on contexts after a full earlier-cohort
regression (1,173,099 contexts). Every preserved coordinate array and all model
summary totals matched. Started paired AA/3Di input preparation with unchanged
profile alignment and masks; retained coverage is not yet established.


Follow-on paired inputs completed and emitted-file checks passed: 88 ready
markers, 84 taxa (73 fungal entries/11 outgroups), 3,986 marker–taxon cells,
843,060 observed paired cells. Full lineage coverage checks passed all 65,750
design cells. The descriptive three-cohort union now covers 525/526 manifest
entries (500 fungal entries/all 25 outgroups), adding 67 newly represented taxa;
Amoeboaphelidium protococcarum F1243177 remains uncovered. This is a union of
separately ready cells, not a merged analysis or proof of 500 unique species.
Same-method cohort integration and predictor-confounding controls remain pending.


Expanded AlphaFold accessibility full audit passed all 13,153 models and
6,910,765 residues (session 22787 exit zero). Complete audit table/model identity,
length, entry-hash and total readback matched the producer and frozen snapshot.
The automatic handoff launched paired-site projection as PID 3397669 under
session 20643; normalization remains queued behind it. This is output
consistency evidence, not independent ASA integration or biological validation.
See metadata/gdm_full_accessibility_audit_receipt.json and
docs/residue-accessibility.md.


AlphaFold accessibility projection completed: 124 markers, 13,066 models and
4,148,852 paired-site observations. All marker observation/missingness counts
matched emitted AA alignments; full projection file hashes passed. Individual
row/ASA readback remains pending. Normalization started as PID 3491107 under
controller session 20643, preserving both reference scales and unclipped values.


AlphaFold paired accessibility normalization completed for all 4,148,852 rows
(controller session 20643 exit zero), preserving both unclipped reference scales.
Started full row-level normalization readback (session 90918), checking source
field/order preservation, inverse scaling, terminal rules and all summary totals.
Producer tables are versioned; independent readback remains running. Projection
against raw residue tables and controlled exposure/evolution tests remain pending.


Full AlphaFold accessibility normalization readback passed every one of
4,148,852 rows and all scale/threshold summaries (session 90918 exit zero).
Raw-residue projection verification and controlled modeling remain pending.

Ecology predictions completed all 675 queued models without interruptions or
memory deferrals (controller session 96562 exit zero). Automatic full artifact
audit launched as PID 3807224 under snapshot-controller session 55711, with
explicit ecology marker links; conversion is queued next. No ecological
association or transition replication has been established by this acquisition.


Ecology artifact audit and conversion completed for all 675 models
(snapshot-controller session 55711 exit zero), with 696 curated marker links
across 18 taxa/81 markers before confidence qualification. Full converted-file
hash checks passed. The full marker panel yields 704 exact-sequence links;
the eight additional links and complete expected grid are retained in versioned
metadata. Launched direct residue mapping (session 37925) after conversion
completed; no queued mapping controller exists. Final link/residue readback and
native confidence qualification remain pending.


Ecology mapping and full residue readback completed: 675 models, 704 exact
marker links across 19 taxa, and all 165,755 matrix-residue links verified.
The expected global links equal all exact full-panel sequence matches and retain
all 696 curated originating links. Started PAE export (session 67992) and native
extraction (session 97473); final confidence/paired coverage remains pending.


Ecology native extraction finished (session 97473 exit zero), using the same
executable, extraction script and native-source hashes as the follow-on cohort.
Full coordinate-feature audit is running (session 25686); initial models pass
reconstruction. PAE export remains live as PID 4052697. Qualification and paired
inputs remain pending.


Ecology PAE export and full exported-file identity/hash readback completed for
675 models (608 MB compressed). The primary coordinate audit completed with
211,633 valid states and 159,854 pre-PAE six-residue confidence-qualified states.
Started PAE qualification after these completed prerequisites. Supplemental
coordinate-array readback subsequently passed all 675 models and 211,633
six-residue confidence minima (session 80243 exit zero). Final qualification
counts and independent PAE-context checks remain pending.


Ecology qualification and full PAE-context readback completed for all 675 models:
211,633 valid contexts checked, 159,639 full-chain residues pass combined feature
pLDDT70/PAE10. Started paired AA/3Di alignment preparation with unchanged profile
matrix, masks and eligibility thresholds. Final usable coverage and ecological
contrasts remain pending.


Ecology paired inputs and emitted-file checks completed: 71 markers, 18 taxa,
645 marker–taxon cells and 134,004 observed paired cells. Full-design lineage
coverage completed. The four-cohort ecological coverage join and all taxon/group
counts matched actual emitted FASTA identities. ESMFold coverage is distributed
across batches: a descriptive same-method union gives 16 common markers across
five Amanita taxa and 26 across three Cenococcum comparison taxa. Proper cohort
integration and evolutionary tests remain pending; these are availability counts.


Combined the three completed ESMFold inventories: all 10,048 sequence/model
identities are disjoint, inference settings match apart from input receipts and
physical GPU identifiers, and all coordinate hashes/provenance passed. Every
emitted inventory object matches its original source exactly. Started full-panel
remapping against the fixed profile matrix; combined qualified encodings and
re-evaluated paired alignments remain pending. The earlier full prediction queue
continues independently; integration uses its frozen audited subset.


Combined ESMFold mapping and qualified encoding integration completed: 10,048
models, 294 mapped taxa, 10,169 marker–taxon links and 2,083,852 matrix-residue
links. Every provenance/mapping row matches the disjoint source union. All
10,048 emitted encoding-summary rows independently match their original source
fields exactly. The unchanged qualified encodings contain 2,009,344 full-chain
residues passing six-residue pLDDT70 and directional PAE10. Started combined
paired alignment preparation on the full 526-taxon, 125-marker design; final
usable coverage and combined branch inference remain pending. Original ESMFold
prediction, follow-on accessibility and full-panel PMSF jobs remain active.


Combined ESMFold paired inputs completed and passed full emitted-FASTA readback:
89 ready markers, 292 taxa (281 fungi and 11 outgroups), 9,332 marker–taxon cells,
22,205 retained marker columns and 1,699,035 observed paired cells. Full-design
coverage still evaluates all 526 working entries and 125 markers (65,750 rows).
The 1,733,521 observed cells in that full coverage table include non-ready or
ineligible cells and are not the emitted-alignment total. Combined ecological
coverage confirms 16 shared markers across five Amanita taxa and 26 across three
Cenococcum taxa; transition replication/review remains unresolved.

Launched 356 combined-cohort fits (89 AA sequence searches, each followed by
three fixed-sequence-topology structural fits). AA searches include 1,000 SH-aLRT
and 1,000 ultrafast bootstrap replicates with bootstrap NNI. Four single-thread
workers with 2 GB per fit; 1–168 hours is a planning scenario, not a measured
runtime bound. Inputs span 4–237 taxa and 73–462 columns per marker. Matched-site
direct geometry for this combined snapshot, uncertainty resampling, reconciliation
and controlled acceleration/coupling analyses remain pending. The full fungal
sampling, larger atlas and all eight evolutionary analyses remain the goal.


Combined PAE manifest integration completed for all 10,048 models. Every
exported matrix and original NPZ checksum matched combined provenance and the
original qualification receipts; every emitted manifest field independently
matched its original source row. The derived manifest preserves original files
and lossless-export/context audits rather than recalculating PAE.

Launched direct geometry for all 89 combined paired markers: up to 704,031 taxon
pairs, on exactly the paired observed sites, with whole-chain/local comparisons
and directional PAE filtering. One CPU worker, 24 GB memory allowance, 10 GB
output headroom; 2–48 hours is a planning scenario. This will benchmark combined
alphabet paths against coordinate changes after supported paired fits complete.
Full numeric/output audit and path comparison remain pending.


Follow-on accessibility calculation completed for all 4,252 models and
1,342,046 residues; its pinned full source audit is running. Queued paired-site
projection and two-reference normalization behind that exact audit controller,
with 17 pinned inputs/dependencies and fresh immutable outputs. The controller
requires a complete matching full audit before proceeding.

Started accessibility for all 675 ecology-cohort models / 212,983 residues, with
the same software, radii, 960 sphere points and 1.4-A probe as the follow-on cohort.
All method configuration fields matched after excluding snapshot and size fields.
Four CPU workers, 8 GiB planning memory and 2 GiB output allowance; .1–3 hours is
a planning range. A pinned full-audit controller is queued. Combined accessibility
integration, independent projected-row verification and controlled surface/core
evolutionary analyses remain pending. Predicted isolated-chain exposure is not
experimental accessibility or evidence of a binding interface.


Added reusable full projected-accessibility readback and started it on the
expanded AF cohort (4,148,852 observed paired cells). It checks each emitted row
against original ASA tables, model provenance, source residue mapping and paired
AA/3Di FASTAs, with exact values and complete observed-cell accounting. This
closes the previously recorded raw-projection verification gap when it passes;
it does not recompute ASA or establish experimental exposure. Execution remains
pending. Follow-on audit/projection and ecology ASA acquisition continue.


Full AF projected-accessibility readback passed: all 4,148,852 cells across 124
markers, 13,066 models and 13,565 marker–taxon cells exactly match raw residue
ASA/confidence/atom/context fields, source mapping and paired alignment identity.
The complete observed grid has no duplicated, omitted or masked emitted cells.
This completes the previously pending raw-projection check, without independently
reintegrating ASA or validating biological exposure.

Follow-on full accessibility audit passed for 4,252 models / 1,342,046 residues;
all emitted audit-model identities, lengths and entry hashes independently match
the source universe. Projection and normalization completed for 843,060 observed
cells in 88 markers / 3,938 models. Tien normalization has zero values above one;
Miller has 2,281 (unclipped). At threshold .25, 36,292 cells (4.3048%) differ
between scale classifications. Full projected-row and normalization readbacks
passed for every one of the 843,060 cells. Ecology projection/normalization is queued behind its exact live
full-audit controller. Controlled evolutionary analyses remain pending.


Completed and independently checked the site-rate/exposure analysis frame for
all 16,571 sites in the earlier 72-marker ESMFold cohort: 132,568 rate values
(four AA/3Di specifications, Gamma and optimized FreeRate), copied exposure
quantiles and copy-review flags, plus amino-acid composition, entropy and
coverage. All source values and fit diagnostics matched; every site composition
and coverage value was independently recomputed from FASTAs. Rates are
conditional phylogenetic estimates, while extant exposure summaries remain
unadjusted. This prepares covariates for controlled coupling/exposure modeling;
it does not constitute an effect test or resolve feature dependence/uncertainty.
The larger combined-cohort supported fits and direct geometry continue running.


Ecology accessibility completed for all 675 models / 212,983 residues. Full
source-coordinate audit passed; emitted audit-model identities, residue counts
and entry hashes match the frozen source universe. Projection and two-reference
normalization completed for 134,004 observed cells in 71 markers / 625 models.
Full projected-row and normalization readbacks passed for all cells.

All three frozen ESMFold acquisition cohorts now have audited accessibility
calculations. Combined raw-accessibility integration and projection remain
pending: the combined alignment has 1,699,035 observed cells versus 1,692,000
in the three separately ready alignment sets, so concatenating their projected
tables would omit newly eligible cells and retain incompatible paired columns.
Reproject original per-model ASA onto the combined alignment after preserving
source audit/configuration provenance. Controlled evolutionary tests remain
pending; the full project scope is unchanged.


Integrated audited ASA for all 10,048 models / 2,738,370 full-chain residues.
All original coordinate/residue checksums and model provenance matched; the
combined audit record explicitly describes a disjoint union of source audits,
not a new coordinate audit. All numerical method settings match; the earlier
cohort used two workers and newer cohorts four, retained as scheduling provenance.
Every emitted audit row and entry hash independently matched the original source
unions, and all 20,096 relative symlinks resolve to the exact source files.

Projection onto the 89-marker combined paired alignment is running. Updated
projection/readback consumers to accept the explicitly named derived-audit
status only with matching union source provenance. The source-audit workflow
continues to use its original status. Normalization, combined row readback and
controlled exposure/coupling analyses remain pending.


Combined accessibility projection and normalization completed: all 1,699,035
observed cells in 89 markers, 9,225 models and 9,332 marker–taxon combinations.
Full projected-row and normalization readbacks passed. Tien has no values above
one; Miller has 5,484, retained unclipped. The .25 threshold differs for 73,826
cells (4.3452%) across scales. These are method sensitivities, not independent
biological transitions or experimentally verified exposure.

Updated predictor-specific coverage union: 22,677 usable marker–taxon cells
across 525 entries (500 fungal entries, all 25 outgroups). Every constituent cell
matched actual emitted AA/3Di FASTA identities. Amoeboaphelidium protococcarum
remains uncovered. This does not establish 500 unique species or authorize
pooling predictor-specific alignments. README now reflects completed ESMFold
cohort integration. Supported combined fits and direct geometry remain active;
controlled exposure/coupling inference and the full project remain incomplete.


Moved the full 513–768-residue prediction queue onto idle GPU0. Fresh exact-sequence
reuse/resource checks retain the same 5,512-sequence source universe, with two
reusable structures and 5,510 prediction candidates across 279 taxa / 108 markers.
The old longer-marker controller was confirmed waiting with no child process,
refreshed input or prediction output, then terminated and marked superseded
before launch. Its terminal receipt prevents accidental restart. Original GPU1
predictions and experimental-control queue continue unchanged.

The replacement controller launched all 5,510 candidates on the explicitly idle
second RTX 6000 Ada, preserving full sequences, four CPU threads and stop-on-OOM
behavior. Updated short-protein observations give 41.92–66.30 inference GPU-hour
scenarios, 99.45 hours including overhead and 37.86 GB quadratic allocated-memory
extrapolation; these are not measured long-protein bounds. Reserve 100 GB output.
No new paid services. Execution, long-protein resource measurements and subsequent
structure/PAE audits remain pending. Reproduce from
`metadata/long_marker_gpu0_prediction_config.json` using
`scripts/run_ready_long_marker_predictions.py`.

GPU0 execution confirmed active alongside GPU1. The first completed full-queue
record is 677 residues with 36.48 seconds inference; this single observation is
not a runtime bound or independent structure-quality validation.


Queued a pinned full audit/conversion handoff behind the active GPU0 longer-marker
producer (PID 2006275, start ticks 181975155). The controller requires exactly
5,510 completed predictions, no OOM deferrals, no interruption and zero remaining
eligible sequences before full PDB/NPZ readback and mmCIF conversion. It validates
matching audited/converted model counts and freezes input/configuration hashes.
Twelve dependencies are pinned. Allowances: one CPU, 8 GB memory, 100 GB output,
.2–12 hours after production; these are planning estimates. An incomplete producer
will stop the handoff for review. Residue mapping and evolutionary integration
remain later stages. Both GPU producers are active.


Prepared the codon-model environment while structural jobs run: HyPhy 2.5.101
from pinned commit 646e16eb4243e5e5588ac340483fcb987552751d, isolated CMake 3.31.6,
GCC 13.3, four build workers and 457 installed-file checksums. The upstream
mitochondrial codon test passed likelihood, omega and total-tree-length assertions
with two CPU threads. Version-specific BUSTED option incompatibilities are
documented from source. No fungal selection or synonymous-saturation results
are claimed; the genus/code screen still needs family-level eligibility work.
See docs/codon-model-environment.md and versioned build/test receipts.


Completed and independently verified within-genus codon diagnostic inputs:
2,084 strict source groups retained in the ledger, 1,712 ready and 372 explicit
coverage exclusions. Every exported codon/AA, taxon membership, 80% occupancy
mask and source coordinate passed full readback (5,728,643 observed taxon-codons).
Ready cases comprise 1,613 standard-code and 99 code-12 groups; 14 ready TFIIB
cases retain copy-reconciliation caveats. Curated-hybrid exclusions affect
14 memberships. Group-specific phylogenetic/divergence inference remains next;
no selection eligibility or positive-selection finding is claimed.


Launched all 1,655 information-screened within-genus nucleotide tree searches
with GTR+F+G4, 1,000 SH-aLRT and 1,000 UFB replicates plus bootstrap NNI;
retain bootstrap trees. Four one-thread/2 GB workers run on the existing host.
The complete 1,712-row information screen reproduces byte for byte; 57 cases
remain explicit information exclusions for this workflow. A timestamped running
observation records 287 completed case receipts. Full tree/report/support audit,
codon divergence and selection eligibility review remain pending. Both ESMFold
GPU queues were confirmed active. See docs/coding-sequence-workflow.md and
metadata/genus_codon_tree_* for configurations, resources and progress provenance.


Independently audited the first 415 completed genus/code nucleotide trees and
all 415,000 saved bootstrap trees. Exact tip grids, branch values, split counts,
model reports, reported total lengths and UFB frequencies passed (rounding
within 0.5 percentage points). Added root-invariance/malformed-tree tests (two
passed). This is an immutable incomplete snapshot, not full-batch completion.
Among 1,686 internal edges, 806 have UFB below 95 and 468 SH-aLRT below 80.
Warning review flags seven cases for saturated nucleotide distances, 18 for
parameter boundaries and 16 for NNI convergence; categories overlap. Retained
127 cases with memory-adjustment warnings separately for investigation. These
observations guide subsequent codon diagnostics, not selection claims.


Investigated the 127,127 repeated IQ-TREE memory warnings against pinned 3.0.1
source. They occur in exactly all 127 four-taxon cases in the audited snapshot.
A slot-count cap of n-2 followed by a minimum int(log2(n)+1) explains the warning
for n=4 despite ample requested memory; reported adjusted allocation ranges
0.011–0.894 MB, and --mem 2G parsing is correct. Archived source files/checksums
and full case-set evidence. No executable or live configuration was changed;
other numerical warnings remain under review. This is a source/log explanation,
not compiled-binary tracing or proof of model adequacy.


Verified all 128 genetic-code entries and both code option names in the running
HyPhy executable, then launched the full 1,655-case global MG94xREV/CF3x4
queue with four one-CPU workers. Cases start when nucleotide-tree execution
receipts become available. The corrected v2 run has 73 completed fits at the
versioned observation; the initial seed-argument failure occurred before model
execution and remains archived. Saved likelihood functions, JSON reports and
logs support independent review. Full tree/fit audits, branch-component scaling,
profile-interval behavior and biological eligibility are still pending. No LRT
or positive-selection result is claimed.


MG94 audit caught and corrected a branch-grid violation: HyPhy's default
zero-branch reduction collapsed internal edges in 48 of 433 completed v2 cases.
Stopped that producer cleanly, preserved all outputs and launched the full
1,655-case v3 queue with explicit --kill-zero-lengths No plus exported-topology
and branch-identity checks. First corrected snapshot: 41 cases/389 branches
passed saved-likelihood reevaluation, codon-grid and numerical-report readback;
maximum likelihood discrepancy 5.46e-12, component additivity error 2e-10.
No numerical flags in this snapshot; full audit and biological eligibility
remain pending. Structure predictions and nucleotide-tree searches continue.


Completed the first full conditional site-coupling grid: 24 regressions across
16,571 sites/72 markers, with marker fixed effects, coverage/confidence controls,
composition/entropy sensitivity, two rate models and two RSA scales. All source
confidence/RSA joins and independent coefficient/cluster-covariance calculations
passed. The sequence-rate coefficient at RSA=.25 is positive in all 24 fits
(23 conditional BH q<.05), but full-model within-marker R-squared is only
5.42–6.90%; no interaction q<.05. A standalone figure, coefficients and methods
are documented in docs/conditional-site-coupling.md. Fixed-tree/rate uncertainty,
marker dependence, extant-exposure weighting and prediction circularity remain
limits; these are exploratory snapshot associations, not final biological tests.


Completed coupling robustness across all 24 specifications: 48,000 paired
marker-bootstrap fits and 1,728 leave-one-marker-out fits, with no singular
bootstrap draws. Full-intercept and absorbed coefficients agree to 6.53e-15;
selected omission/bootstrap fits independently match expanded-row least squares.
Observed-mean sequence-rate contrasts remain positive and RSA contrasts negative
under every single-marker omission. All unadjusted marker-bootstrap intervals
for those contrasts exclude zero; interaction intervals all include zero.
These are conditional resampling sensitivities, not new multiplicity-adjusted
biological tests. Added observed-mean contrast figure and full results in
docs/conditional-site-coupling.md. Full 1,655 nucleotide trees finished; their
full bootstrap audit is running.


All 1,655 nucleotide trees and corrected MG94 fits are now complete. Full tree
readback passed 1,655,000 bootstrap trees and 6,721 internal edges. Full MG94
readback passed 1,655 saved likelihood reevaluations and 18,407 branch records,
with no numerical review flags (maximum LL discrepancy 9.10e-12; component
additivity error 2e-10). Separate source binding verifies all expected case IDs
and exact audited input topologies. Full warning review retains 33 saturated-
distance, 53 boundary, 83 NNI-convergence and two high-gap cases; categories
overlap. Memory warnings match exactly all 505 quartet cases. Codon-specific
saturation, biological eligibility and uncertainty propagation remain pending.


Combined ESMFold direct geometry and readback are complete: all 704,031
within-marker pair records checked (686,583 accepted, 17,448 excluded), plus
89 independently recomputed geometry/PAE examples selected by hash. The 89-marker
paired branch fits remain running, so fitted-path benchmarking remains pending.
Also pinned upstream FitMG94's equilibrium-opportunity normalization formula
for dS/dN; this source review does not yet produce normalized divergence tables.


Completed site-opportunity normalization for all 1,655 MG94 fits and 18,407
branches. Independent checks first detected an upstream stop-compaction defect:
GTT/TAC/TCT S/NS entries were overwritten in both codes (12 scalar errors).
An isolated one-line correction passed all 244 opportunity assertions and exact
patch-reproduction checks; installed libraries and fitted models remain intact.
Every case's frequency-weighted opportunities independently agree between Python
and HyPhy (maximum discrepancy 2.23e-15), and all normalized branch outputs passed
inverse readback. Three regression tests pass. Very large dS estimates occur
(maximum 99.09); identifiability/saturation and alignment review remain essential
before selection inference. Normalization uses the documented equal-alternative
opportunity convention, not a new branch-specific omega model.


Completed longest-dS-branch conditional likelihood slices for every genus/code
case: 11,585 evaluations and 1,655 restored-baseline checks. No sampled branch
multiplier improves the saved fit beyond numerical precision (9.10e-12).
Largest estimates can still have peaked fixed-nuisance curves; this does not
establish absence of synonymous saturation. Added a full case review linking
target branches to taxon splits, codon coverage and existing copy/tree warnings,
plus a figure of the 12 largest estimates. Nuisance-reoptimized profiles and
alignment/model eligibility reviews remain pending.

Full raw-log/table readback of the longest-branch slices passed for all 1,655 cases and 11,585 points. The saved figure and scripts accompany this checkpoint. Structural fit queues remain active: 420/496 AF-based and 68/356 combined ESMFold fits at this observation; these counts are execution receipts, not final audited results. Both ESMFold GPU queues and the species-tree/orthology jobs remain running.

### Expanded tree-distance/geometry handoffs queued

Queued checksum-bound controllers for both running expanded paired-fit datasets:
124 AlphaFold-covered markers and 89 combined ESMFold markers. Successful full
fit completion triggers the full fit audit, direct-geometry/tree-path benchmark,
descriptive rank summary and figure generation. Both real-source preflights and
three prerequisite-failure tests passed; both controller processes were verified
live and waiting on the exact original fit producers. This advances automation
of the expanded structural-branch benchmark; downstream results remain pending.
See `docs/tree-path-geometry.md` and the versioned controller configurations.

### Nuisance-reoptimized branch-parameter profiles running

Launched the full 1,655-case MG94 diagnostic grid: seven fixed values of the
longest branch parameter plus one unconstrained reoptimization per case, totaling
13,240 fits and fresh saved-model likelihood readbacks. A first-launch exporter
failure was caught by readback, diagnosed and retired with all artifacts retained.
The corrected v2 queue now checks the complete branch parameter grid; initial
completed cases reproduce their optimized likelihoods in fresh processes.
The queue is live. Full summary/audit, interpretation and eligibility review
remain pending; these parameter profiles are not dS intervals or selection tests.

### Full selected-assembly contamination-report inventory and audit completed

All 519 selected NCBI assembly-version FCS reports were retrieved. Validation
retains 472 publisher-MD5-verified reports and 46 reports verified by identical
repeat retrieval but lacking publisher MD5 entries. One actual publisher-checksum
disagreement persists on repeat (F104408/GCA_013423385.1), and seven external
genomes remain outside this NCBI report assessment. Among the 518 usable reports,
86 assemblies have EXCLUDE/FIX/TRIM regions; 4,039 action-preserving region rows
are retained. Full original raw-row readback and independent interval-union
checks passed. No genomes or proteins were removed. CDS/structural-marker overlap
and biological review are the next dependencies; these observations do not
establish contamination-free assemblies or confirmed contamination in every
protein. See `docs/assembly-quality-workflow.md` and `metadata/selected_assembly_fcs_*`.

### Contamination-report coding overlaps and analysis exposure audited

All 4,039 validated FCS regions were intersected with 4,479,554 CDS feature rows
across 86 assemblies and independently audited. The 28,150 intersections involve
11,092 annotated proteins (11,062 retained representatives). Fifty-four original
marker observations overlap flagged coding regions: 44 N. cerealis, nine
A. colombiana and one review-only A. queenslandica observation. Exact identity and
input-presence checks locate 16 exposures in the earlier ESMFold paired data,
19 in combined ESMFold, none of these 54 in expanded AlphaFold paired data, and
29 affected N. cerealis genus-codon cases. Source-quality review and explicit
exclusion sensitivities now take priority over interpreting these affected inputs
as biological acceleration. Baseline results are retained; no running input or
inference script was changed. Details and reproducible commands are in
`docs/assembly-quality-workflow.md`.

### FCS-overlap sensitivity refits launched for all affected ESMFold markers

Prepared and independently read back both complete affected-marker sensitivity
sets. Omitted 16 observations across 16 earlier ESMFold markers and 19 observations
across 18 combined ESMFold markers under the EXCLUDE/FIX/TRIM CDS-overlap rule.
All affected markers remain eligible, no columns become all-missing, and every
retained character matches the baseline. The 56/71 unchanged ready-marker
complements keep their original inputs and estimates. Launched 64 and 72 new
supported sequence/structural fits respectively; both producer processes were
verified live. Full fit audits, topology/path comparisons and downstream
site-rate/coupling sensitivities remain pending. Baseline runs are preserved.

### Full MAFFT guide audited and remaining crossed PMSF runs queued

The 526-taxon, 63,750-site MAFFT LG+F+G4 guide completed and passed full readback.
The two guides share 470/523 internal splits; their normalized unrooted RF
distance is 0.10134. Full independent edge-split checks passed; path association
is descriptive only. Composition concerns persist (526/526 nominal failures for
MAFFT). Both guides are unsupported and not final species phylogenies.
A live, pinned controller now queues the other three alignment/guide PMSF
combinations serially after successful completion of the current first run,
retaining existing resource checks and the exclusive PMSF lock.

Full FCS-omission codon accounting also completed: all 29 affected cases would
have three remaining taxa and fall below the existing four-taxon gate; the other
1,626 are unchanged by this specific omission. The affected cases are not
refitted under the current design and are not interpreted as null results.

### FCS sensitivity direct-geometry subsets audited; fit handoffs queued

Prepared exact baseline-geometry subsets for every remaining pair in the changed
sensitivity markers: 108,005 earlier-ESMFold pairs and 280,954 combined-ESMFold
pairs including exclusions. Every output field was checked against baseline
measurements with only the cohort taxon count updated. Fresh full grid/mask/
character/eligibility audits passed, plus 16/18 independent sampled numerical
geometry checks. Both fit-to-benchmark controllers were verified live and waiting
on the exact sensitivity fit producers. Full fit audits and sensitivity path
summaries will follow successful fitting; matched-pair baseline comparisons,
uncertainty and updated site-rate coupling remain pending.

### Alternative-hit recovery inventory and exploratory gene-tree review

Reviewed all raw HMMER outputs and final BUSCO calls for the 54 FCS-overlapping
marker observations: 1,004 marker–protein hits, with no unflagged alternative in
final BUSCO calls. Raw outputs yield 41 alternatives across 35 observations that
lack recorded coding overlap and reach the absolute dataset score cutoff. Their
exact sequences and audit metadata are staged for review; 37 have exact CDS
translations and four non-triplet CDS records. They are not accepted replacements.
The installed BUSCO relative-score filter was documented without attributing an
unverified rejection cause to individual candidates.

For one exploratory affected marker (129234at2759), the completed gene tree
contains a split grouping N. cerealis with four Candida entries, reported SH-aLRT
87.2. This marker-specific discordance contrasts with concatenated-guide
neighbors and motivates source review; it is not contamination confirmation or
whole-genome taxonomic reassignment. Baseline/sensitivity inputs remain unchanged.


Alternative candidate alignment checkpoint: 41 proteins across 33 markers were
projected onto the unchanged profile-matrix columns using the original pinned
HMMs. Independent Stockholm/coordinate readback checked 14,798 candidate site
rows, 9,855 observed candidate residues and 5,757,596 unchanged baseline
characters. Retained-site coverage min/median/max is 19.97%/74.63%/100%.
All candidates have observed residues; 1,465 all-missing baseline marker rows
are explicitly identified for exclusion before review-tree inference. None of
these candidates is an accepted replacement. See assembly-quality-workflow.md
and metadata/fcs_alternative_profile_alignment_{receipt,audit}.json.
Original and long-sequence ESMFold producers, both baseline paired-fit queues,
both FCS sensitivity fit queues, first PMSF run and nuisance-reoptimized codon
branch profiles were verified live at this checkpoint. The project remains
incomplete; gene-copy placement and the queued sensitivity results are pending.

Alternative-copy phylogenetic review advanced: prepared all 33 affected marker
alignments under the existing per-sequence gene-tree coverage rule. 38/41
candidates qualify; three remain below threshold and explicitly recorded.
Independent input readback checked 17,399 source rows and 5,284,204 retained
characters, including exact retained/dropped identifier sets. Started 33
exploratory copy trees with two concurrent IQ-TREE jobs (two threads, 4 GB each,
MFP LG/WAG/JTT+F+Gamma, 1,000 SH-aLRT, identical tips preserved). Producer PID
3765266/session 65513 and both initial IQ-TREE children were verified live;
revalidate before any restart. Planning estimate about 20.8 hours from completed
baseline trees, range 8–72 hours; no new paid infrastructure. Inputs/configuration
and source/model/software hashes are pinned. No candidate orthology or
replacement has been accepted; full tree audits and copy placement remain pending.

At this checkpoint, original expanded AlphaFold paired fits had 460/496 individual
fit receipts, combined ESMFold 120/356, and FCS sensitivity queues 20/64 and 8/72.
All four producer processes were verified live and no full root receipt existed.
The baseline marker-tree queue had 48 completed individual trees. These are
partial completion counts, not fully audited final analyses.

FCS omission exposure covariates completed: exact retained ASA projections for
16 earlier ESMFold and 18 combined ESMFold markers contain 223,145 and 427,510
observed residues. Full independent raw-source readbacks passed for every row,
including protein/model/residue identity, confidence, atom count, context,
AA/3Di character and observed-cell coverage. Both reference normalizations and
all summary totals passed separate full readbacks. At relative ASA 0.25,
normalization-scale disagreements are 4.14% and 4.07%, respectively, descriptive
convention sensitivity only. Full command manifest and receipts committed under
metadata/fcs_accessibility_execution.json and metadata/*fcs_accessibility*.
Rate refits and exposure aggregation across retained sensitivity taxa are still
needed before updated sequence–structure coupling tests. Earlier and combined
cohorts overlap; this checkpoint is not independent replication.

Matched-site FCS exposure comparison completed for all affected markers:
earlier ESMFold 16 markers/2,492 sites/2,167 omitted residue observations;
combined ESMFold 18 markers/2,950 sites/2,667 omitted observations at 2,502 sites.
Exact retained-row equality and complete observed-taxon alignment grids passed.
All exposure/confidence quartiles were checked by independent linear
interpolation and NumPy. Largest absolute median-RSA changes are 0.03954/0.03967
(Tien/Miller) earlier and 0.01856/0.04695 combined. Extant covariate shifts are
descriptive and do not establish robustness of the pending rate/coupling refits.
Receipts and full marker summaries: metadata/*fcs_site_exposure*.

Queued complete changed-marker rate/exposure continuations for both FCS
sensitivity cohorts. Each eleven-stage controller waits for its exact live
fit-audit/geometry parent, verifies completion and source pins, then performs
Gamma4 and FreeRate exports/audits, all-fit four-start optimization checks,
selected-fit comparison/readback, retained-taxon parsimony/exposure summaries
and analysis-frame construction. 384 likelihood fits are planned for earlier
ESMFold (16 markers), 432 for combined (18 markers). Existing resources only,
four single-thread workers/2 GB each per inference stage; fresh 64 GiB memory
and 20 GiB disk headroom gates. Both preflights and live wait states verified.
Controllers PID 3786355/session 6692 and PID 3786372/session 76128; revalidate
before any restart. Full configurations, commands, source/dependency hashes
and resource estimates are versioned under metadata/fcs_rate_analysis_*.
Merging unchanged baseline markers and re-estimating coupling remain pending.

Marker-tree support snapshot v4 now covers 48/125 completed trees and 22,705
internal edges (22,454 reported SH-aLRT values, 251 unreported). An independent
source/input/split readback reconstructed all retained sequences from the frozen
profile matrix, verified coverage exclusions and 8,465,734 alignment characters,
and recovered all internal splits through undirected graph-edge removal.
Fitted model labels, attached support/branch values and report hashes checked.
77 markers remain pending; completion-order bias prevents full-batch discordance
claims from this snapshot. See docs/phylogenetic-workflow.md and versioned v4
snapshot/readback receipts and marker tables. Baseline FCS flags remain relevant.

Full species-tree FCS sensitivity inputs and guide runs advanced. Both v2
matrices retain all 526 taxa (including all 25 outgroups), 125 markers and the
original retained column sets. 53 observations across 48 markers in two fungi
are masked (44 Naganishia cerealis, nine Acaulospora colombiana). Full readback
checked 25,788,202/33,532,500 profile/MAFFT characters. All taxa retain observed
residues. A v1 affected-taxa receipt-labeling defect was corrected before any
inference; matrix bytes were unchanged and retirement evidence is retained.
Launched two homogeneous LG+F+G4 guide sensitivities on existing resources with
16 threads/32 GB each. Parent processes verified live; launch identities and
resource estimates in metadata/fcs_species_guide_*. Full tree audits, topology
comparisons and heterogeneous-model sensitivities remain pending.

Dating evidence advanced: verified the existing publisher-MD5 timetree archive,
inventoried all files, and audited four annotated mean chronograms against all
node-age/interval summaries. Each contains 153 tips/305 nodes; 1,216 mean-age
edge differences agree within 4.55e-13 Ma. Their 182-entry name maps contain 29
entries absent from each tree, now explicitly marked. 57 actual published tips
per scenario have exact-name project matches requiring identity review. No
dates were transferred: these summary trees do not supply joint posterior
branch-duration samples. Added calibration qualification and uncertainty
requirements in docs/dating-workflow.md, grounded in primary timetree and
calibration-provenance literature. Dating and per-time structural estimates
remain pending, without blocking relative-divergence analyses.

Dating calibration inventory advanced: downloaded and hashed the primary
supplementary calibration document, catalogued all 27 entries across 24 numbered
groups (21 published minimum bounds, six soft maxima), and retained separate A/B
entries. Independent PyMuPDF extraction checked every heading's page, numeric
age and bound type against the pdftotext-derived catalogue. Four entries have
both anchor names exactly matched to project names; other candidates need
explicit clade mapping, not automatic rejection. No fossil calibration was
approved or date transferred. Full provenance and qualification limits are in
docs/dating-workflow.md and metadata/published_calibration_*.
Alternative-copy tree producer was verified live; no completed copy trees were
available at this checkpoint.


Repository status reconciliation: updated the README and leading milestone table
against current receipts and live job identities. Corrected stale MAFFT-guide,
FCS-inventory and accessibility states; exposed the FCS caveat beside the initial
coupling result; linked dating evidence and retained all eight objectives as
incomplete. Historical execution entries remain unchanged and date-specific.

Dating input implementation review completed for pinned McmcDate source commit
cce22b458fbb583863ca53b0f44793719363980c. Verified that CSV node targeting uses
MRCA tip names and provided age boundaries require probability-mass fields.
The literature heading catalogue alone is not an executable calibration model.
Recorded the inspected duplicate-checking scope and need for cross-node
compatibility review. The exact study software version and executable fungal
calibration table remain unresolved; no priors were silently filled in and no
clock run was launched. Source files and tree inventory are archived outside
Git; URLs, hashes and findings are in metadata/dating_input_source_review.json.

Alternative assembly availability review completed for the two taxa with
EXCLUDE-overlapping marker observations. Both frozen fungal GenBank/RefSeq
catalogs and fresh, fully paginated NCBI Datasets taxon queries return only the
selected accession for each species: GCA_030039265.1 (Naganishia cerealis) and
GCA_910592055.1 (Acaulospora colombiana). Saved-response readback verified all
artifact hashes, exact taxids and accessions. No alternative accession was found
within this scope; external archives, historical records and mislabelled taxa
are not ruled out. No assembly replacement or sequence-provenance resolution
was claimed. Reproducible inventory and evidence are in
scripts/inventory_fcs_alternative_assemblies.py and
metadata/fcs_alternative_assembly_inventory_*.json. Alternative-copy trees and
FCS omission sensitivities remain in progress. At this checkpoint, individual
paired-fit receipts reached 488/496 expanded AlphaFold, 160/356 combined ESMFold,
60/64 earlier ESMFold FCS sensitivity and 28/72 combined FCS sensitivity;
complete-batch audits remain pending. Conditional MG94 branch-parameter
profiles reached 907/1,655 complete cases; the producer remains live.

Full FCS-sensitivity normalized exposure inputs completed for both ESMFold
cohorts. Reused all 56/71 unchanged markers and verified the affected 16/18
markers against their audited normalized rows. Full retained alignment grids
contain 712,769/1,696,368 observations, with 2,167/2,667 omissions. Independent
readback checked every retained field against the baseline after the specified
whole marker/taxon omissions (11,404,304/27,141,888 values); source and output
hashes all match. No site-rate frame or coupling results were merged yet.
Expanded AlphaFold paired fitting advanced to 492/496 individual receipts;
its parent and downstream controller remain live, without a complete-batch
receipt. The previous turn was progress (new scoped assembly-availability
evidence and reproducible inventory); this turn completed the full exposure
inputs needed by the forthcoming coupling sensitivity.

Two supported paired-fit batches completed and passed their full audits:
expanded AlphaFold 124 markers/496 fits/26,758 paired branches and earlier
ESMFold FCS omission 16 markers/64 fits/3,428 paired branches. Near-zero branches
occur in 491/496 and 63/64 fits, respectively; no branch is at least 10. These
point estimates therefore require uncertainty analysis, not branch-length
ratios. Their completed direct-geometry benchmarks contain 737,851 and 105,834
accepted pairs, respectively, with 4,960/640 independent path checks. Both
three-panel rank figures were visually inspected; receipts and visual reviews
are now versioned in metadata. Earlier FCS G4 site rates and their audit have
completed; its controller has advanced to R4, with optimization and coupling
refits still pending.

Compared earlier FCS paths against baseline fits on exactly the same 105,834
retained accepted pairs: 423,336 path values and all 20 fixed fields per pair.
All 448 reported before/after rank correlations independently reproduce through
SciPy within 1.11e-16. The authoritative comparison is
results/phylogeny/paired-path-fcs-comparison-esmfold-v2; v1 was superseded after
the independent check detected default decimal-parser effects on geometry ties.
Explicit round-trip parsing resolves the discrepancy without changing trees or
input geometry. Median baseline/sensitivity path-rank agreement is 0.988 for AA
and 0.968–0.982 for the structural models, but individual marker/model agreement
can fall to 0.638. Geometry correlations on matched pairs also vary, with changes
roughly -0.127 to +0.085. These descriptive checks do not establish robust
coupling or acceleration. Scripts and compact marker-level tables are versioned
as compare_fcs_paired_paths.py, readback_fcs_path_ranks.py and
metadata/esmfold_fcs_path_comparison_*.

The first of 33 alternative-copy trees completed: marker 4764044at2759,
476 tips/391 columns. Audited artifacts and exact tip/input grid, then verified
all 952 candidate/selected focal-to-tip paths through independent graph
traversal (maximum difference 1.33e-15). The original N. cerealis protein lies
on a five-tip split with four Candida species; its alternative copy lies on a
three-tip split with N. friedmannii and N. liquefaciens. Both separating edges
report SH-aLRT 100; this is not a contamination probability or acceptance of
orthology. The remaining queue remains live. Full findings, source-quality
caveats and next qualification requirements are in docs/fcs-copy-tree-review.md
and metadata/fcs_copy_tree_review_4764044at2759.json. No protein replacement or
new prediction was made from this single-tree result.

Launched conditional paired resampling for all 124 expanded AlphaFold markers
after their completed fit audit. The established workflow uses 200 paired draws
at circular block lengths 1, 10 and 30: 74,400 draws and up to 148,800 fixed-
topology fits (AA LG+F+G4 and structural AF+G4). Prelaunch resource plan reserves
32 GB memory and 250 GB disk with a broad 2–96 hour planning envelope; available
host resources were checked and no paid resources provisioned. Producer PID
3921215/session 20841 and its start identity are recorded in metadata. Startup
inspection found eight active fit children and 141 completed draw receipts.
The full audit controller is live under session 38481 and will wait for the
exact producer before checking the complete output. Script/source hashes,
configuration, resources and launch identities are in
metadata/gdm_expanded_paired_resampling_*.json. No resampling batch is declared
complete. These draws measure conditional sampling sensitivity, not calibrated
confidence intervals: topology, dating, alignment, prediction and model
uncertainty remain unpropagated, and local blocks do not preserve all nonlocal
structural-alphabet feature dependencies. The preceding turn was a verified
wait with advancing live jobs; this turn started the next required uncertainty
analysis on the completed AlphaFold dataset.

Earlier ESMFold FCS workflow completed all 256 optimization refits, full audits,
selected Gamma/FreeRate comparison, full rate readback, exposure and changed-
marker frame. Merged its 16 changed markers with 56 unchanged baseline markers:
all 16,571 sites/132,568 rates retained, 712,769 observed taxon-site exposures.
Reran all 24 coupling specifications and independently verified every coefficient
and cluster covariance (maximum errors 1.03e-14/4.03e-16). All sequence-rate
coefficients remain positive; 20/24 meet the same 72-test BH threshold versus
23/24 at baseline. No matched focal coefficient changes sign; all interaction
intervals include zero. Threshold changes are not coefficient-difference tests.
The full 72-pair coefficient comparison and reviewed figure are versioned.
Rewrote the current sensitivity section in docs/conditional-site-coupling.md;
FCS marker-resampling/influence checks and broader uncertainty/circularity
controls remain incomplete, as does the combined-cohort sensitivity workflow.

FCS coupling marker-resampling/influence checks completed: 48,000 bootstrap
fits and 1,728 single-marker omissions, no singular bootstrap fits, all point
and selected expanded-row solves verified. Full marker order and 2,000 sampling
multiplicity vectors equal baseline; all 120 saved percentile intervals checked.
Observed-mean sequence-rate contrasts stay positive and observed-mean RSA
contrasts negative under every single-marker omission, with unadjusted bootstrap
intervals excluding zero. All interaction intervals include zero, and one basic
AF/FreeRate/Tien interaction changes sign under a marker omission. These are
conditional sensitivity results, not new BH tests or resolved phylogenetic/rate
uncertainty. Reviewed figure and full numerical/provenance records are versioned.

Second alternative-copy tree reviewed: 345792at2759, 484 tips/226 columns,
two alternative copies and 1,936 independent focal path checks. Selected
A. colombiana copy lies with Serendipita entries; its alternative lies near
Diversispora/Ambispora. The N. cerealis selected/alternative contrast recurs,
with weak support for some finer relationships explicitly retained. No candidate
replacement or contamination conclusion accepted. Details and support caveats
are in docs/fcs-copy-tree-review.md and the marker-specific review receipt.

Full MG94 nuisance-profile execution completed: 1,655 cases, 13,240 optimized
points and fresh saved-fit readbacks. New full audit checked 66,200 artifact
hashes, 226,696 fitted parameter values, all unchanged model text, likelihood
replays and the exact aggregate grid. Maximum likelihood readback discrepancy
8.64e-11; additional readback of 26,480 target-t/omega summary fields agreed
exactly. Initial audit-parser attempts incorrectly included the fixed reference
exchangeability among optimized variables; corrected its scope and verified
that fixed model content remains identical. No fitted outputs were changed.
Four Malassezia cases have a constrained grid likelihood above the unconstrained
reoptimization (maximum 1.995), so restart diagnostics remain necessary.
All 29 FCS-exposed cases are flagged. Full results and caveats are recorded in
docs/coding-sequence-workflow.md and metadata/genus_mg94_branch_parameter_profile_*.
This is an audited finite-grid diagnostic, not dS intervals, global-optimum
proof, saturation clearance or selection inference.

Launched 32 unconstrained restarts for the four flagged MG94 cases, using each
case's eight saved profile solutions. Each completed restart releases only
the target constraint, preserves the remaining model definition and undergoes
a fresh saved-likelihood readback. Two cases have completed; the full run is
still active. Source profiles remain unchanged. Resource estimates and the
reproducible runner are recorded in metadata/genus_mg94_flagged_restart_plan.json
and scripts/restart_flagged_genus_profile_optima.py. Finite restarts do not
establish a global optimum; revised profiles remain a subsequent task.

September 16 completion update: all 32 restarts and fresh likelihood replays
have now finished. A subsequent complete readback verified 194 artifact hashes,
model exports, saved likelihoods and summary arithmetic. All four cases found
better likelihoods than the previous best evaluated solution (gains 0.0112 to
2.2559), with large target-t estimates. Revised profiles and identifiability
assessment remain outstanding; no selection eligibility is inferred. Results
are archived in metadata/genus_mg94_flagged_restart_* and documented in
docs/coding-sequence-workflow.md.

Combined ESMFold FCS omission fits and geometry benchmark completed for all
18 affected markers: 72 fits, 6,058 paired branches, 275,252 accepted pairs,
5,702 exclusions and 720 independent graph path checks. All four downstream
receipt hashes and 11 listed artifacts were independently rechecked; the
three-panel figure was visually reviewed. All 72 fits contain near-zero
branches, so branch ratios are inappropriate. Descriptive within-marker rank
correlations are not acceleration tests or adjusted sequence–structure
coupling. Site-rate optimization is running through the existing controller;
the full combined-cohort coupling analysis remains pending. Receipts are
versioned under metadata/esmfold_combined_fcs_sensitivity_*.

Reviewed a third completed alternative-copy tree, 4803900at2759: 419 tips,
106 columns and 838 independent paths. The selected N. cerealis copy groups
with three Candida entries, while the candidate groups with N. liquefaciens.
Short-alignment, long-branch and source-quality caveats are retained in
docs/fcs-copy-tree-review.md; no candidate replacement was accepted.

User-requested temporary cooling applied September 14: current project process
trees share 24 CPU cores; the active ESMFold process alternates 30-second
run/rest periods. No hardware cap was available without administrator access.
Initial live checks found all 913 examined thread affinities restricted and
GPU power falling to 83–89 W during rest; predictions continue. This exception
expires after 48 hours and does not change future defaults. Recovery, provenance
and caveats: docs/temporary-cooling-20260914.md. Earlier GPU ETA is superseded.

September 14, approximately 09:10 Eastern: user requested 3.5 hours at prior
compute levels. Stopped the cooling controller cleanly and verified original
CPU affinities across 896 live threads, with GPU prediction running normally.
A one-off user systemd timer is active and waiting to restore the same cooling
settings at approximately 12:40 Eastern today. Original cooling expiry remains
September 16 00:56 Eastern. Schedule and checks are versioned in
metadata/temporary_full_compute_20260914_schedule.json; this is not a permanent
change to project defaults.

DGX Spark feasibility assessment completed September 14. Existing SSH access,
10-Gb/s direct link, idle GB10/20 ARM cores/~116 GiB available memory confirmed.
SLURM GPU smoke test and four ESMFold checks completed on Spark (job 20912,
5m06s, exit 0). Checks used unchanged runner/weights in an isolated project
venv; inherited PyTorch 2.10 differs from the workstation environment. Full
output readback passed, but one 700-residue prediction differs by 6.417 A
whole-chain CA RMSD (1.646 A on jointly confident positions). Spark predictions
therefore remain separate from the atlas pending source-compatibility review.
Four tested predictions were 4.1–4.9 times slower per sequence on Spark.
No production work moved; local cooling schedule unchanged. CPU offload remains
an option using native ARM tools. Details: docs/dgx-spark-resources.md.

September 16 reboot recovery: the user-requested pause suspended 32 processes,
but the scheduled continuation has no success record. Host rebooted at 08:58
Eastern; all prior process identities and the transient timer were gone.
Recovered eight reviewed restartable jobs, retaining 2,477 saved ESMFold
prediction receipts for startup validation and reusing checked tree/search
outputs. Four incomplete domain tables were archived before rerunning their
chunks. PMSF had no serialized checkpoint: preserved old logs, checked all
input pins, and began the same analysis in a fresh recovery directory.
The 24-core and 30-second GPU run/rest limits remain active. OrthoFinder
continuation and replacement completion controllers are still pending review.
Details and recovery provenance: docs/recovery-20260916.md.


September 16, 10:40 Eastern: verified the temporary compute increase remains
active at 48 physical cores and 90-second GPU run / 10-second rest cycles.
GPU0 is computing at its default/max 300 W limit, around 80 C without reported
thermal slowdown; GPU1 is idle. Prediction, PMSF and isolated OG0000017 FastTree
processes are live. Both site-rate wrappers and downstream controllers remain
live; their terminal completion is not yet asserted.

Matched experimental controls now have reviewed protein-balanced tables and a
figure: 39/36/16 proteins eligible at joint pLDDT cutoffs 0/70/90. Independent
aggregation verified 125,515 metric values and 216 quantiles. All-protein median
paired ESMFold-minus-AlphaFold experimental RMSD differences are +0.6154,
+0.4688 and +0.2995 Å. The same-16-protein comparison explicitly retains the
caveat that residues and entries can change. These are descriptive selected-set
agreements, not training-independent accuracy estimates. The 35 longer control
sequences remain pending. Details: docs/experimental-structure-workflow.md;
figure: docs/figures/experimental_predictor_controls.svg.


September 16, 10:43 Eastern: GPU1 now runs the next 12 experimental-reference
controls (513–741 residues), disjoint from every existing prediction queue.
Both GPUs were observed at 100% utilization during processing. The new service
pins the unchanged prediction/audit/conversion scripts and advances only after
all 12 predictions complete. Twenty-three longer controls remain deferred for
memory planning; no truncation or completed-control claim. Plan and launch:
metadata/experimental_control_long_tier_{plan,launch}.json.


September 16, 10:46 Eastern: scheduled 20 additional full experimental controls
(788–1,063 residues) behind the 12-control GPU1 tier. The live handoff requires
successful predecessor prediction/audit/conversion before checking GPU idleness
and launching. Memory planning uses observed 768-residue peaks with explicit
extrapolation/headroom; OOM requires review. Three original controls remain
above the new ceiling (1,361/1,468/2,413 residues). Documentation and resource
provenance: metadata/experimental_control_next_tier_*. No completion claim for
the pending tier.


September 16: archived the completed 124-marker AlphaFold Gamma-rate audit
(496 fits, 178,876 rate rows). Combined diagnostics cover 852 Gamma export fits
across separate source cohorts. One ESMFold AA export, marker 5005750at2759,
is 2.5645 log-likelihood units below its original fit and requires targeted
Gamma optimization review; queued FreeRate refits do not resolve that baseline.
Original outputs remain unchanged. See docs/site-specific-evolutionary-rates.md.


September 16: two stricter Gamma refits of ESMFold marker 5005750at2759 completed.
The better observed fit improves the original by 5.7239 log units; full tip,
edge and site-grid readback passed. A revised immutable Gamma baseline and
affected-marker sensitivity remain necessary before final coupling claims.
The original live controller inputs remain intact. Also verified completion
and conversion of all 12 longer experimental controls; the next GPU1 batch
started automatically. Receipts are versioned in metadata.


September 16: revised ESMFold Gamma baseline passed the full 356-fit/88,820-rate
audit. Byte comparisons confirm only marker 5005750at2759 AA outputs changed.
Replaced the verified idle downstream controller before any stage started;
its successor uses this audited baseline for the existing planned 1,424
FreeRate optimization diagnostics and exposure-frame workflow. The FreeRate
producer remains active. No duplicate diagnostic execution or final coupling
claim. See docs/site-specific-evolutionary-rates.md.


September 16: all 12 newly predicted longer experimental controls now have
completed three-way structural comparisons and independent readback: 759
threshold rows, 387 accepted, 372 excluded and 1,161 superpositions checked.
All 12 proteins qualify without confidence filtering, nine at joint pLDDT70,
and one at90. The unfiltered median paired ESMFold-minus-AlphaFold experimental
RMSD is +5.5277 Å. Small selected cohorts, domain orientation and experimental
context prevent general accuracy or biological interpretation. Full tables
and limitations: docs/experimental-structure-workflow.md.


September 16: completed all annotated Pfam-region comparisons for the 12 longer
experimental controls: 35 hits, 2,220 threshold rows, 1,368 eligible. Independent
grid/count/eligibility readback passed; dual-superposition checks passed. In the
exploratory Q04305 case, individually fitted annotated regions agree locally
(ESMFold ~0.4–0.7 Å) despite ~19–26 Å discrepancies under whole-mask alignment.
This supports treating whole-protein placement separately from within-region
change; it does not establish biological domain motion or functional novelty.
Context/PAE review remains pending. Details: docs/experimental-structure-workflow.md.


September 16: completed region-level ESMFold PAE summaries across 12 longer
controls, with 237 rows and 584 independently checked quantiles. In the Q04305
example, within-region median PAE is ~0.9–1.4 Å versus ~24–28 Å between WD40
repeats and UTP15_C, consistent with uncertain relative placement despite high
focal confidence. PAE is predicted uncertainty, not measured error or biological
motion. Also archived the completed 89-marker/356-fit ESMFold FreeRate audit;
the revised downstream controller has started its 1,424 optimization fits.


September 16: scheduled the 1,361- and 1,468-residue experimental controls behind
the current 20-control batch, using measured GPU allocation through 1,017
residues and an explicit allocator-reservation caveat. The live handoff requires
completed predecessor prediction/audit/conversion and a fresh idle GPU.
Sequences remain full length and disjoint from all other queues. One 2,413-
residue control remains deferred for a memory-saving execution approach;
its estimated allocation exceeds 48 GiB. No completed-results claim for queued
controls. Plans and evidence: metadata/experimental_control_large_tier_*.


September 16: queued a separate attention-chunk16 ESMFold attempt for the last
2,413-residue experimental control, plus an intentional 741-residue method
comparison. Local code review confirms the supported attention setting but
not a guaranteed total-memory reduction; an OOM stops for review. The live
handoff waits for the two larger controls to complete validation first.
Original predictors/settings remain available; variant results require an
explicit numerical comparison before integration. Source review and resource
plan: metadata/experimental_control_chunk16_*. No completed prediction claimed.


September 16: added a reviewed, reproducible Q04305 region-placement figure
showing full-sequence PAE alongside matched regional/whole-mask RMSD for both
predictors. All 16 plotted values match the archived table; rendering and file
hashes checked. Added the executed matched-control and region/PAE procedures
to the methods draft, retaining exploratory-selection and biological-context
limitations. Figure: docs/figures/control_region_placement.svg.


September 16: completed marker-to-full-guide branch projection for both source
cohorts: 45,155 marker edges, 90,310 guide-specific rows. All 426 independent
pruned trees and 67,994 compatible path sums passed readback. Internal edges
mapping uniquely to the same full split in both guides number 1,975/9,065 for
ESMFold and 4,600/13,193 for AlphaFold. Collapsed paths and discordance are
explicit; no structural change is assigned to an arbitrary full-tree branch.
Supported-tree variants and final acceleration models remain pending. Details:
docs/marker-species-edge-projection.md. Also verified all 20 longer GPU controls
passed artifact audit/conversion; the next two larger controls are running.

### September 16: remove temporary compute throttling

At the user’s request, terminated the temporary resource controller cleanly. Its restoration receipt reports no errors: GPU 0 no longer receives scheduled 10-second pauses, and tracked CPU processes/threads regained their original affinity. Both GPUs were observed at 100% utilization with normal 300 W power limits after restoration. Existing job-specific thread counts and memory limits remain; hardware thermal protections and scientific prediction settings are unchanged. This applies to the current execution, not a global resource policy. See `metadata/compute_unthrottled_20260916_receipt.json`.

### September 16: 77-control experimental comparison snapshot

Completed the 20-protein, 788–1063-residue control batch: 968 accepted and 595 excluded chain/model/threshold rows. The coordinate readback independently checked all 1,563 residue masks and 2,904 superpositions. A reusable pandas summary audit checked 19,825 hierarchical metric values, 216 cohort quantiles and all 60 disposition rows. Eligible protein counts are 18/16/1 at joint predicted confidence cutoffs 0/70/90. Median paired RMSD differences (ESM–experiment minus AF–experiment, chain-level differences aggregated hierarchically) are +3.147194/+1.741811/−1.475913 Å; the final value represents only one protein.

Combined the disjoint 45-, 12-, and 20-protein batches into `results/experimental_structures/predictor-controls-summary-combined77-v1`. Source tiers and all 231 protein/threshold dispositions remain explicit. Independent readback checked exact source-table preservation and all 864 combined/tier cohort quantiles, counts and delta signs. Across the 77-protein universe, 69/61/18 qualify at cutoffs 0/70/90, with median paired RMSD differences +1.208679/+0.646558/+0.299482 Å. Restricting to the same 18 proteins gives +0.543479/+0.471100/+0.299482 Å. These are descriptive results for selected experimental references, not general accuracy estimates, independent biological replicates or causal confidence/length effects.

The 1,361- and 1,468-residue controls also completed prediction, artifact audit and conversion; their geometry comparisons remain pending. Thus 79 of the original 80 control targets have completed standard-setting predictions. The separate chunk-16 experiment completed its intentional 741-residue repeat, but the 2,413-residue target recorded `oom_deferred`; the controller stopped with `failed_requires_review`. No automatic retry or truncation was attempted. The successful repeat still needs artifact validation and a direct comparison against its chunk-64 counterpart before interpretation. GPU 0 production and the phylogenetic analyses continue; GPU 1 became available after this bounded attempt. Reassigning its next disjoint production batch is the next scheduling task.

### September 16: second GPU reassigned to longer marker production

Launched 4,363 full-length 769–1,024-residue marker predictions on GPU 1 in `fungal-extended-marker-gpu1-20260916.service`. Candidates link to 276 taxa and 99 markers; the full selected band includes 277 taxa because six exact-sequence models were reusable. Preparation froze the current AFDB retrieval snapshot and checked all existing local prediction directories, including the experimental controls. Independent prelaunch readback verified all candidate sequence hashes, length eligibility, source artifact hashes, complete marker links and zero overlap with GPU 0’s 513–768-residue queue. No sequence is truncated.

The resource plan uses 20 measured long-control timings: quadratic/cubic scenarios predict 111–138 GPU-hours, with a conservative 207-hour planning allowance including 50% overhead. This is not a deadline or confidence interval. Planned resources are one exclusive 48 GiB GPU, four CPU threads, 64 GiB host memory and 100 GiB output allowance; launch required 200 GiB free disk. Both GPUs use continuous processing with normal hardware power/thermal protections. The staged runner stops for review on OOM or incomplete predictions; after successful prediction it runs full artifact readback and mmCIF conversion automatically. The transient service does not survive reboot; its launch PID/create time and input/command hashes are recorded in `metadata/extended_marker_gpu1_*`. Species-wide structural coverage and downstream evolutionary estimates still require integration after completion.

The preserved 741-residue chunk-16 repeat passed full artifact readback (`experimental-controls-chunk16-partial-audit-v1`); its failed 2,413-residue companion remains one eligible unfinished target. The audit flag `partial_prediction_snapshot: false` denotes terminal-chunk rather than live-snapshot mode, not completion of both targets. Direct chunk-16 versus chunk-64 numerical comparison remains pending.

### September 16: attention-chunk numerical control

The audited 741-residue repeat has exactly identical serialized Cα coordinates and native confidence/PAE arrays at attention chunk sizes 64 and 16. Independent NumPy Kabsch and Bio.SVDSuperimposer calculations agree; RMSD is numerical zero across all residues and joint native confidence cutoffs 70/90 (741/494/223 residues). No residues cross confidence cutoffs 50/70/90. Configuration comparison verified that checkpoint, model source, packages, precision, recycles, seed and GPU agree; only the reviewed runner, queue, eligibility length limit and chunk size differ. Audited prediction-table membership and receipt/artifact hashes bind both inputs.

Peak allocated GPU memory was 12,297,761,280 versus 12,297,509,376 bytes: only 251,904 bytes (0.24 MiB) lower at chunk 16. Inference times were 52.90 versus 48.79 seconds, an uncontrolled single-run comparison. This does not establish a general speed advantage, useful memory reduction or numerical equivalence on other proteins. Source inspection confirms the changed setting feeds triangle-attention chunking; it does not chunk every tensor operation. The 2,413-residue companion still failed with OOM. No more retries are justified by this control alone. The duplicate 741-residue result is a method check and is not added as another protein to the benchmark.

The two completed 1,361/1,468-residue controls have exact AlphaFold/ESMFold crosswalks across 17 experimental entries. Their matched geometry, independent coordinate readback, hierarchical summaries and independent summary check are running sequentially under the pinned `experimental_control_large_geometry_plan.json`. GPU 1's extended marker batch has now emitted its first verified prediction (1,003 residues; 114.71 seconds).

The two-control comparison subsequently completed: all 51 masks and 81 superpositions passed independent coordinate readback, followed by the independent hierarchical summary check. Both proteins qualify at cutoffs 0/70 and neither at 90. The new disjoint 79-protein snapshot (`predictor-controls-summary-combined79-v1`) has 71/63/18 eligible proteins at cutoffs 0/70/90. Independent readback verified eight source tables, 237 dispositions and 1,080 cohort quantiles/counts/signs. Four length tiers and the unchanged 18-protein common cohort remain explicit. The original 77-protein snapshot is preserved; the failed 2,413-residue target and duplicate chunk-setting control are excluded from the 79 distinct proteins.

### September 16: numerical phylogeny–structure branch integration

Assembled and independently read back `projected-gamma-branch-frame-v1`: 45,155 marker edges, 180,620 paired numerical branch estimates and all 90,310 prior guide mappings retained. The corrected ESMFold Gamma baseline and expanded AlphaFold Gamma fits remain source-separated. All 852 tree edge-length sets were reconstructed by graph traversal; all unique-full-edge marker coverage groups were checked. Unique internal edges number 188 ESMFold and 139 AlphaFold, with 71/81 respectively having at least ten markers. No collapsed-path change was allocated to individual species branches. These coverage and conditional branch estimates prepare the branch-level model; they do not establish acceleration, positive selection, time rates or a supported species phylogeny. Full numerical frame stays outside Git; source pins, readback, scripts and coverage summaries are versioned.

Revalidated live production: both GPU prediction PIDs, revised ESM optimization controller (four worker children plus wrapper), AlphaFold rate exporter, PMSF inference and OG0000017 FastTree repair remain live. The family repair log advanced to 38,780/55,596 top-hit sequences, so no restart is warranted.

### September 16: supported marker-tree and guide conflict integration

Advanced the audited supported marker-tree snapshot from 48 to 70 completed trees. Independent input/split readback checked 14,170,535 alignment characters and 33,068 internal edges. The remaining 55 trees remain pending. Added descriptive SH-aLRT-cutoff (80/95) agreement/conflict diagnostics against both homogeneous guides, retaining coverage failures, unresolved cases and explicit strongest-conflict witnesses. The 146,440 rows include repeated restricted splits from collapsed guide paths and cannot be treated as independent evidence or full-batch concordance factors. Full matrices stay outside Git; compact counts, source pins and reusable scripts are versioned. Supported species topology, complete-marker discordance and reconciliation remain required.

### September 16: expanded coupling execution handoffs

Installed and verified live waiting controllers for the revised 89-marker ESMFold and 124-marker AlphaFold frames. After the ongoing rate-model audits complete, they will automatically execute 24 coupling fits per cohort, 48,000 marker-bootstrap fits per cohort and all 2,136/2,976 leave-one-marker-out fits. Source identities, counts, resources, scripts and dynamic plan construction are pinned; failure or an incomplete upstream handoff stops for review. No coupling results are claimed yet. Both GPUs were observed at 100% utilization. See the expanded-handoff section of docs/conditional-site-coupling.md for resource estimates and remaining inference limitations.

### September 16: domain-profile detection mapped onto marker trees

Completed two-policy minimum-change mapping for detectable Domain/Family/Repeat Pfam profiles on 70 audited marker trees. Preserved all 2,620 marker/profile/policy dispositions; 185 baseline and 81 conservative combinations qualified (five detected/five undetected tips). Recorded 6,631 possible-change edges while retaining all optimal-reconstruction ambiguity. Exhaustive small-tree enumeration verified the numerical algorithm; independent empirical readback verified 126,667 tip states and all 266 minimum scores. Among 81 paired eligible combinations, 46 scores fall under overlap/partial-hit masking. This quantifies annotation sensitivity rather than establishing biological domain gains/losses, rooted lineage events, copy-number changes or full architecture evolution. Full-proteome annotation and reconciliation remain pending. Source pins, methods, resource plan and profile summaries are versioned.

### September 16: combined experimental benchmark figure

Completed and visually reviewed the 79-control, six-panel benchmark figure, retaining four length tiers and confidence-filter exclusions. Independently checked 152 plotted protein/threshold rows, 456 RMSD values, 12 coverage counts, six medians and common-protein membership. Corrected draft legend overlap and used shared scales before publishing the SVG. Updated README marker-tree and paired-fit status from the latest receipts. Full prediction, supported-phylogeny and coupling workflows remain active; this figure does not complete the broader objective.

### September 16: additional full-proteome hit annotation in production

Pinned all 52 completed Pfam chunks (6,498,460 hits) and launched streaming annotation across the additional 5,654,720 query sequences. The original four HMMER workers remain live on subsequent chunks. At the checkpoint, 18 annotation shards totaling 2,472,034 hits had completed; the full snapshot remains in progress. An independently implemented raw-field readback is queued on the exact annotation producer identity. Resource plans, launch identities and parsers are versioned. No full-proteome architecture, overlap resolution or domain-absence claims are implied.

### September 16: 6.5-million-hit annotation snapshot independently checked

Completed streaming annotation and independent raw-field readback for all 52 pinned chunks: 6,498,460 hits. Every row retained source identity/order, coordinates, scores, E-values, posterior accuracy and HMM coverage. Pfam descriptive metadata remains from the pinned parser and was not independently reparsed. The remaining 12 search chunks are still pending. Installed and verified a live completion service to annotate/read back those chunks after successful raw-search completion and assemble a disjoint full64 catalog. Synthetic catalog gates reject missing, duplicate and modified shards. This completes the current annotation snapshot, not the entire full-proteome annotation, marker-partition merge, architecture or evolutionary analysis.

### September 16: fungal–outgroup separation and root uncertainty

Verified that both 526-tip guides contain the 501-fungal/25-outgroup separating split, without rooting either saved tree. Independent graph traversal checked both guides and all 70 audited marker trees. Only 16 marker trees contain the exact corresponding split; at SH-aLRT cutoffs 80/95, supported agreement/conflict/unresolved counts are 14/50/6 and 9/36/25. All marker split and support values cross-match earlier audited tables, and both guide root partitions yield identical diagnostics. This documents material discordance for subsequent reconciliation and ancestral work, not a finalized root or a biological explanation. Complete-marker and supported-mixture sensitivities remain necessary. Revalidated both GPU producers, full Pfam workers, rate controllers, PMSF and family-tree repair as live during this turn.

### September 16: expanded AlphaFold FreeRate milestone

Completed and audited all 496 FreeRate fits across 124 markers (178,876 site-rate rows). Independently read back the initial Gamma/FreeRate comparison: 178,876 site and 107,032 branch comparisons. Five lower-likelihood FreeRate cases remain explicit. The existing controller advanced automatically to the 1,984-fit optimization stage; four IQ-TREE workers were verified live. Initial summaries and flags are archived separately so optimized results can be compared against them. Expanded coupling remains queued behind optimization and frame validation.

### September 16: independent full-proteome query linkage check launched

Launched `readback_full_domain_protein_links.py` to independently reconstruct all 5,815,847 representative protein links across 526 taxa from source FASTA files. The readback checks sequence hashes, exact taxon/protein/query-source identity, disjoint search partitions, complete additional-query coverage and all per-taxon counts. It uses a parser separate from input preparation. The single-CPU resource plan and live PID/create-time/command are recorded under `metadata/full_domain_protein_link_readback_*`; output is `results/domains/full-protein-link-readback-v1`. The process remains live and has not yet produced a completion receipt. This is a prerequisite for preserving species-specific proteins when joining domain annotations, not a completed architecture analysis. Both GPU producers, rate controllers, full Pfam workers, PMSF and family-tree repair were independently verified live; GPU utilization remained 100% on both cards.

### September 16: full-proteome query linkage independently verified

Completed the source-FASTA readback for all 5,815,847 protein links and all 526 taxa in 453 seconds. Every accession, taxon, sequence identity and query partition matched; all 5,654,720 additional sequences were represented. The 58,879 reused marker sequences correspond to 60,033 protein links; four marker queries are outside the representative-proteome universe. Receipt hashes and all taxon/aggregate totals were checked before archiving. Synthetic re-pinned fixtures rejected missing, duplicate, wrong-taxon and wrong-partition links. This completes the linkage prerequisite for domain annotation joins, not full architecture evolution. Live downstream controllers and their script pins were revalidated; rate optimization and Pfam search remain in progress.

### September 16: automatic validated family-tree installation queued

Added and launched a pinned handoff behind the exact live `OG0000017` rebuild wrapper. It requires successful isolated completion, verifies unchanged inputs, independently checks source tips and explicit finite nonnegative branches with OrthoFinder's bundled ETE-derived parser, then atomically replaces only the audited empty production tree. Six disposable fixtures verified successful installation and rejection/preservation behavior for invalid or incomplete inputs. The live service is waiting; no production tree has been replaced yet. Resource estimates, scripts, native parser hashes and launch identities are versioned. Full native from-trees continuation, reconciliation and biological validation remain separate outstanding stages.

### September 16: native orthology continuation prerequisites resolved further

Native source review and executable file-locator readback found that the interrupted full-run log lacks required cluster/tree metadata and both species-tree artifacts. Prepared a separate explicitly constructed continuation descriptor; the native reader now resolves the audited clusters and existing tree directory without altering the original log. The descriptor remains not launch-ready: IDs species-tree preparation and isolated file-write behavior require review. Installed defaults also enable MSA/tree output updates and a second ortholog-inference call, so a new results directory alone is insufficient evidence of source protection. Archived the reproducible preparation script, native before/after readback and exact source hashes. No premature reconciliation launch occurred; live prediction, tree repair and other analyses continue.

### September 16: native reconciliation restart executed in disposable fixtures

Executed actual installed OrthoFinder continuation on two four-taxon/three-family software fixtures. A supplied species-name tree without the native IDs tree causes an internal missing-file failure but returns exit code zero. Preparing the corresponding rooted IDs tree allows `--from-trees --no-fix-files` to complete without changes to the source fixture files. Independent readback verified all 13 genes, three HOGs/resolved trees, 42 directed ortholog pairs, one expected terminal duplication and the species root split. Archived reproducible fixture generation, artifact validation, resource plan, commands and source/output hashes. This resolves a concrete restart dependency and establishes that full-run completion requires output validation beyond process return code. Full-cohort species-tree staging and reconciliation remain pending.

### September 16: both full-cohort guide sensitivities prepared for reconciliation

Prepared rooted taxon-name and native-ID inputs for both audited 526-tip homogeneous guides. The 501-fungal/25-outgroup split is explicit, with equal root-edge halves used solely as a serialization convention. Native conversion checked IDs/topology; high-precision intended files preserve original unrooted branch lengths rather than the converter's rounded diagnostic output. Independent graph deletion checked all six files, full bifurcation, exact IDs/root split and 6,294 canonical-edge comparisons. These are conditional guide alternatives, not final supported phylogenies. Original trees remain unchanged; isolated staging, completed repair, reconciliation and supported-guide/root sensitivity analyses remain pending.
