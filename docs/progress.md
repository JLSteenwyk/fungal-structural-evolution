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
