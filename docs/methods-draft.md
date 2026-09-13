# Methods draft: completed data preparation

This draft describes work executed so far. Phylogenetic inference, structural clustering and evolutionary tests are not yet completed and are not represented as results here.

## Sampling and sequence acquisition

The working analysis manifest contains 501 fungal species and 25 non-fungal outgroups. Candidate selection used annotated public assemblies and a deterministic taxonomic diversity design, supplemented by published proteomes from underrepresented lineages. The current fungal circumscription follows the broad NCBI classification; sensitivity to a narrower fungal circumscription remains planned. A further fungal candidate, Saccharomyces jurei, lacks an available protein dataset and is retained in the candidate manifest with an exclusion record. Species selection, ecological curation and quality exceptions remain subject to review.

We acquired 519 proteomes from assembly-specific NCBI files and seven from published Figshare bundles. Publisher checksums, local SHA256 digests, source versions and retrieval locations were recorded. Original files were preserved. QC inputs removed terminal stop/period markers and excluded internally interrupted or noncanonical sequences according to the recorded normalization policy. The published Abeoforma proteome contains terminal periods, which were tracked separately before normalization. QC inputs retain alternative products; they are not gene-collapsed proteomes.

## Broad marker recovery

BUSCO 6.1.0 was run in protein mode with the pinned eukaryota_odb12.2 dataset (125 markers; dataset creation 2026-05-13), using four threads per job and at most four simultaneous jobs. All 526 jobs completed successfully. Per-taxon summaries were checked for internal count consistency. These scores describe recovery of a broad eukaryotic marker panel and do not establish genome completeness independently of evolutionary marker loss, divergence and annotation quality. No uniform exclusion threshold was imposed. Microsporidia have a median complete-marker recovery of 42.4%; Pirum gemmata and Abeoforma whisleri each recover 5.6%. These observations motivate lineage-specific QC and sensitivity analyses, not automatic removal.

Single-copy complete hits were extracted only from successful jobs whose input checksums matched the recorded proteomes. Marker sequences were checked for exact protein accession and sequence equality against source inputs. The full extraction contains 59,840 sequences across 125 markers and all 526 screened taxa, with missing and duplicated markers retained as explicit absences from the single-copy matrix. MAFFT 7.525 alignment was launched using `--auto --inputorder`, four concurrent alignments and two threads per alignment. Validation requires identical taxon IDs and ungapped residues. Alignment filtering, phylogenetic model choice and tree inference remain pending.

## Annotation and isoform reconciliation

All 519 assembly-matched NCBI GFF downloads passed checksum and feature-format validation. Protein accessions were followed through annotation Parent links to gene features. Multipart CDS rows were treated as one feature, while multiple gene associations and absent parent links were retained as unresolved. Four published outgroup GFFs use transcript IDs matching protein IDs; the Creolimax GTF uses transcript identifiers and, for 50 additional products, exact gene identifiers. Phaffia rhodozyma's deposited annotation has CDS/mRNA features but no gene features; its products therefore remain unresolved for gene-copy inference.

For the two sanchytrids, published protein-header coordinates were checked against deposited genome contigs using exact translation. All 7,220 Amoeboradix proteins and 9,367 of 9,368 Sanchytrium proteins were verified. Contig-name aliases were accepted only when a unique node/length prefix and exact translation agreed. The remaining Sanchytrium interval exceeded the deposited contig length. Verified ORFs are provisional loci and are not treated as resolved gene/isoform models.

A reproducible representative baseline chooses the longest protein per uniquely mapped gene, breaking ties lexically by accession. Alternative products remain in source proteomes with selection audit records. Unresolved and provisional ORF products are retained independently and flagged against direct use as gene-copy counts. Full-dataset representative preparation is in progress; sensitivity to representative selection remains required.

## Structural acquisition and remaining analyses

Existing structure candidates are nominated by exact full-sequence matches to UniProt entries with AlphaFoldDB cross-references. Current AFDB metadata and coordinates are retrieved with version and source provenance. Accepted models require exact agreement between input, API and CIF polymer sequences and complete alpha-carbon residue coverage. Per-model confidence summaries are retained. Retrieval is ongoing; downloaded models have not yet passed the additional domain, PAE, orthology and prediction-source assessments needed for evolutionary interpretation.

Species-tree inference, gene-tree reconciliation, domain annotation, missing-structure prediction, structural comparisons, all eight evolutionary analyses, and mechanistic case studies remain outstanding. No claim of structural acceleration, adaptation, ecological association or functional novelty follows from the current data-preparation results.

### Executed exploratory direct comparisons

After the preparation described above, a structural mapping snapshot linked 452 marker proteins across 118 taxa to 159,837 matrix positions. Direct paired comparisons used matched unambiguous amino acids and Cα coordinates at pLDDT thresholds 50, 70 and 90 in both models. Pairs required at least 50 qualified residues and coverage of at least half their shared profile positions. Proper-rotation superposition RMSD and symmetric local Cα distance changes were computed, alongside uncorrected sequence differences at the same sites; full definitions and provenance are in the structure-integration workflow.

The executed batch yielded 865 qualifying within-marker taxon pairs across 111 markers at one or more thresholds, including 858 at pLDDT ≥70. These are exploratory geometry measurements, not phylogenetically adjusted tests or branch rates. No significance test was performed. Domain/PAE checks, source effects, prediction circularity and the eight evolutionary analyses remain outstanding.
