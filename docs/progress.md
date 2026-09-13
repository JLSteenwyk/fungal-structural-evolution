# Progress and evidence

| Milestone | State | Completion evidence required |
|---|---|---|
| Local repository and original objective | Created | Git history and docs/objective.txt |
| GitHub remote | Published and public visibility verified via GitHub API | https://github.com/JLSteenwyk/fungal-structural-evolution |
| Catalog discovery | NCBI fungal catalogs inventoried; additional sources pending | metadata/source_receipts.json and metadata/catalog_summary.json |
| 500 fungi + 25 outgroups | Not selected | Reviewed manifest with coverage audit |
| Data QC and structure inventory | Pending | Per-taxon reports and file checksums |
| Species tree and discordance | Pending | Alignments, trees, support and sensitivity |
| Families, domains and reconciliation | Pending | Memberships, trees and reconciliation |
| Structural atlas | Pending | Structures, provenance, clustering and annotation |
| Evolutionary objectives 1–8 | Pending | Estimates, uncertainty, tests and sensitivity |
| Figures, methods and case studies | Pending | Reproducible artifacts and validation proposals |

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
