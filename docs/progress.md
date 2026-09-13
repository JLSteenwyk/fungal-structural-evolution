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
