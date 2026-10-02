# Methods draft: data preparation and exploratory structural comparisons

The October 2 [expanded covariance handoff](full-expanded-covariance-20261002.md)
computes complete new species factors and endpoint/entity incidence for every
expanded logical case. It preserves distinct genes sharing versioned models,
repeated targets/controls and original selection membership. The full design
has 9,812 exact zero-sum taxon patterns on 526 working-kernel tips. Independent
shared-entity graph and original tree-edge reconstruction gate acceptance.
The full handoff passed independent readback and provenance closure:
1,052,632 entity occurrences, 3,330 family components, rank 301 and all
481,376,720 ordered pattern-pair covariances across five working trees. Closure
checked 1,830,546 bindings and both original journals. No new variance
components or calibrated evolutionary effects have been fitted.

The October 2 [expanded case/measurement handoff](full-expanded-matched-measurements-20261002.md)
adds complete source-linked logical case identities and all target/control
order comparisons for the expanded fixed matching. The full case producer and
independent SQL reader passed 4,250,692 selections and all unmatched partitions;
75,188 logical cases retain gene/family/taxon/gene-tree-node and physical model
identities separately. The directed catalog passed all 1,123,936 states through full provenance
closure. The case/mask join passed independent Decimal readback and full
provenance closure for all 150,376 rows, with 2,369,783 source/artifact bindings
and both original journals. Both-orders-only means
and complete four-order envelopes preserve missing and excluded states for
RMSD and explicitly native TM dissimilarity. Separate Decimal reconstruction
checks derived arithmetic. These are dependent measurement records for future
expanded covariance and model fitting, not calibrated effects or completed
phylogenetic adjustment. See the workflow for all source, mask, denominator,
quarantine, resource and original-journal requirements.

This draft describes executed methods through 22 September 2026. Early exploratory datasets remain explicitly identified by their snapshot sizes. Full structural-atlas construction and the final phylogenetically integrated evolutionary analyses remain unfinished; completed preparation and conditional estimates are not treated as final biological results.

## Sampling and sequence acquisition

The working analysis manifest contains 501 fungal taxon entries and 25 non-fungal outgroups. A name-based identity audit flags 21 incompletely identified fungal labels. Subsequent assembly-linked literature curation identifies two hybrids, including one not detected by the original name-only screen; these 501 entries are not yet established as 501 distinct accepted species. Candidate selection used annotated public assemblies and a deterministic taxonomic diversity design, supplemented by published proteomes from underrepresented lineages. The current fungal circumscription follows the broad NCBI classification; sensitivity to a narrower fungal circumscription remains planned. A further fungal candidate, Saccharomyces jurei, lacks an available protein dataset and is retained in the candidate manifest with an exclusion record. Species selection, ecological curation and quality exceptions remain subject to review.

We acquired 519 proteomes from assembly-specific NCBI files and seven from published Figshare bundles. Publisher checksums, local SHA256 digests, source versions and retrieval locations were recorded. Original files were preserved. QC inputs removed terminal stop/period markers and excluded internally interrupted or noncanonical sequences according to the recorded normalization policy. The published Abeoforma proteome contains terminal periods, which were tracked separately before normalization. QC inputs retain alternative products; they are not gene-collapsed proteomes.

## Broad marker recovery

BUSCO 6.1.0 was run in protein mode with the pinned eukaryota_odb12.2 dataset (125 markers; dataset creation 2026-05-13), using four threads per job and at most four simultaneous jobs. All 526 jobs completed successfully. Per-taxon summaries were checked for internal count consistency. These scores describe recovery of a broad eukaryotic marker panel and do not establish genome completeness independently of evolutionary marker loss, divergence and annotation quality. No uniform exclusion threshold was imposed. Microsporidia have a median complete-marker recovery of 42.4%; Pirum gemmata and Abeoforma whisleri each recover 5.6%. These observations motivate lineage-specific QC and sensitivity analyses, not automatic removal.

Single-copy complete hits were extracted only from successful jobs whose input checksums matched the recorded proteomes. Marker sequences were checked for exact protein accession and sequence equality against source inputs. The full extraction contains 59,840 sequences across 125 markers and all 526 screened taxa, with missing and duplicated markers retained as explicit absences from the single-copy matrix. MAFFT 7.525 alignment was launched using `--auto --inputorder`, four concurrent alignments and two threads per alignment. Validation requires identical taxon IDs and ungapped residues. All 125 MAFFT alignments completed and passed identity and residue-preservation checks. They contain 697,866 full-protein alignment columns; requiring unambiguous residues in at least half the taxa present in each marker retained 63,750 columns in a 526-taxon alternative matrix. Separately, HMMER 3.4 profile alignments against the pinned BUSCO models completed for all 125 markers, yielding 50,941 profile match columns. Requiring unambiguous residues in at least half the taxa present for each marker retained 49,027 columns in the 526-taxon concatenation. Both full-taxon homogeneous guide trees and all 125 supported marker trees have completed inference and validation. Two of four crossed C20-PMSF species-tree executions have passed full saved-profile/tree/bootstrap readback. Their current scope and the remaining species-tree requirements are described below; the homogeneous guides are not final species-tree results.

## Annotation and isoform reconciliation

All 519 assembly-matched NCBI GFF downloads passed checksum and feature-format validation. Protein accessions were followed through annotation Parent links to gene features. Multipart CDS rows were treated as one feature, while multiple gene associations and absent parent links were retained as unresolved. Four published outgroup GFFs use transcript IDs matching protein IDs; the Creolimax GTF uses transcript identifiers and, for 50 additional products, exact gene identifiers. Phaffia rhodozyma's deposited annotation has CDS/mRNA features but no gene features; its products therefore remain unresolved for gene-copy inference.

For the two sanchytrids, published protein-header coordinates were checked against deposited genome contigs using exact translation. All 7,220 Amoeboradix proteins and 9,367 of 9,368 Sanchytrium proteins were verified. Contig-name aliases were accepted only when a unique node/length prefix and exact translation agreed. The remaining Sanchytrium interval exceeded the deposited contig length. Verified ORFs are provisional loci and are not treated as resolved gene/isoform models.

A reproducible representative baseline chooses the longest protein per uniquely mapped gene, breaking ties lexically by accession. Alternative products remain in source proteomes with selection audit records. Unresolved and provisional ORF products are retained independently and flagged against direct use as gene-copy counts. Representative preparation completed for all 526 taxa, retaining 5,815,847 proteins and recording 111,898 alternative products. Sensitivity to representative selection remains required. OrthoFinder 3.1.5 inference began with a diversity-selected 64-taxon computational core (54 fungi and 10 outgroups), followed by full assignment and guide-defined clade discovery. Both expanded family partitions now include all 5,815,847 proteins exactly once; their annotation linkage has passed independent validation. Reconciliation under alternative guides is running and final reconciled orthology remains unvalidated.

## Structural acquisition and remaining analyses

Existing structure candidates are nominated by exact full-sequence matches to UniProt entries with AlphaFoldDB cross-references. Current AFDB metadata and coordinates are retrieved with version and source provenance. Accepted models require exact agreement between input, API and CIF polymer sequences and complete alpha-carbon residue coverage. Per-model confidence summaries are retained. Retrieval is ongoing; downloaded models have not yet passed the additional domain, PAE, orthology and prediction-source assessments needed for evolutionary interpretation.

Full-proteome domain annotation and five local prediction cohorts are complete. Final species-tree selection, validated gene-tree reconciliation, remaining structural acquisition, full-atlas structural comparisons, completion of all eight evolutionary analyses, and validated mechanistic case studies remain outstanding. No claim of structural acceleration, adaptation, ecological association or functional novelty follows from the current data-preparation results.

### Executed exploratory direct comparisons

After the preparation described above, a structural mapping snapshot linked 452 marker proteins across 118 taxa to 159,837 matrix positions. Direct paired comparisons used matched unambiguous amino acids and Cα coordinates at pLDDT thresholds 50, 70 and 90 in both models. Pairs required at least 50 qualified residues and coverage of at least half their shared profile positions. Proper-rotation superposition RMSD and symmetric local Cα distance changes were computed, alongside uncorrected sequence differences at the same sites; full definitions and provenance are in the structure-integration workflow.

The executed batch yielded 865 qualifying within-marker taxon pairs across 111 markers at one or more thresholds, including 858 at pLDDT ≥70. These are exploratory geometry measurements, not phylogenetically adjusted tests or branch rates. No significance test was performed. Domain/PAE checks, source effects, prediction circularity and the eight evolutionary analyses remain outstanding.

### Taxonomic identity review

A deterministic screen of all 526 labels flagged 21 names containing `sp.`, `cf.` or `aff.`, and one explicit hybrid, Saccharomyces cerevisiae x Saccharomyces kudriavzevii. The review table preserves assembly accessions and flags pending species identity assessment and exclusion sensitivities. No taxa were removed from the current production runs. Absence of a label flag does not verify accepted taxonomy or rule out hybrid ancestry. Final species-level sampling and bifurcating species-tree interpretation require this review.

### Executed PAE sensitivity

Version-matched PAE matrices were retrieved and validated for all 422 distinct models in structural mapping snapshot-v2. For the 858 direct comparisons qualifying at pLDDT ≥70, the original matched residue sets, superposition RMSDs and local-neighborhood pair counts were independently reproduced. Confidence at each PAE threshold (5, 10 and 15 Å) required both directional PAE entries in both models to pass. Mean absolute Cα distance changes and residue-pair counts were calculated for confident and uncertain subsets of all nonadjacent pairs and of local pairs within 15 Å in either model; pairs must be separated by at least three sequence positions in both proteins. No PAE data were missing. At the 10 Å threshold, the median fraction passing among all nonadjacent pairs was 0.7455; the median across comparisons of their mean absolute distance change was 0.4849 Å unfiltered and 0.3064 Å among confident pairs. Different pair composition precludes treating this difference as an effect estimate. The complete manual review ordering is descriptive and does not define independent events or statistical discoveries. PAE is predicted alignment uncertainty and does not provide calibrated error bars on these distances.

### Executed alignment correspondence comparison

Profile and MAFFT alignments were compared using identical protein-residue identities retained by both methods' 50%-occupancy masks. Cross-taxon residues assigned to the same column define a residue-pair edge. Counts of profile edges, MAFFT edges and their intersection were obtained from a contingency table of column memberships; edge Jaccard is intersection divided by union. The calculation does not enumerate all pairs. Complete-audit status, alignment hashes, source-protein agreement, Stockholm residue coordinates and recomputed occupancy rules were verified. The common retained universe contains 22,269,691 residues; pooled edge Jaccard is 0.925532 and the median across 125 markers is 0.944903 (range 0.573155–0.993988). Per-column coverage and agreement are retained for later sensitivity analyses. Shared agreement cannot establish alignment accuracy, and these results do not yet assess topology or evolutionary-rate robustness.

### Domain annotation search design

Pfam 38.2 profiles and associated entry-type, clan and active-site resources were retrieved from a version-specific EBI archive using pinned publisher MD5 values and local SHA256 provenance. Active-site transfer has not been implemented. All 30,134 HMMs passed accession, length and gathering-threshold checks. Full marker sequences were reconciled against the complete 526-taxon extraction; 59,840 protein records collapsed to 58,883 exact unique sequences for search, preserving every taxon/protein association. All profiles were partitioned into 64 deterministic chunks balanced by summed HMM match-state lengths. HMMER 3.4 searches used the whole unique marker-sequence target database with curated sequence and domain gathering thresholds (`--cut_ga`), four concurrent jobs and two HMMER worker threads per job. Marker and additional full-proteome searches subsequently completed; the merged annotation and candidate-architecture procedures are described below. The prepared parser retains raw overlaps and Pfam entry types rather than equating each hit with a discrete structural domain or resolved architecture.

### Assembly contiguity and deposited metadata

Assembly-version-matched NCBI statistics reports were acquired for all 519 NCBI inputs and validated against publisher checksums and accession headers. Whole-assembly and primary-assembly metrics were kept separate; missing fields were not imputed. Whole-assembly contig N50 was available for 474 reports, while 45 fungal reports omitted it. Sequencing method, deposited assembly type/representation, reported coverage and source BioSample/BioProject links were retained. All seven external genomes were measured directly from their checksum-verified deposited FASTAs using record lengths, base composition and ambiguity counts. FASTA-record N50/L50 were computed without gap splitting or organelle filtering and were not treated as equivalent to NCBI contig metrics. All 526 taxa were joined to validated broad BUSCO summaries; contiguity and recovery were plotted descriptively in separate NCBI/external panels. These measurements do not establish contamination status, biological ploidy or annotation correctness.

### Ecological candidate evidence and initial species curation

The publisher-hosted FungalTraits genus supplement was checksum-pinned and matched by exact genus labels to the selected fungal entries. It supplied candidate records for 482 entries; 19 lacked exact matches and the 25 outgroups were left unassigned. Source primary and secondary lifestyles, other trait fields, spreadsheet coordinates and taxonomy were preserved without promoting genus-level information to confirmed species traits. Six species-level statements from primary genomic studies were recorded separately, with source locators, assembly identity and taxonomic/strain scope limitations. These include two Amanita species whose published asymbiotic classifications conflict with a genus-level ectomycorrhizal projection. Confirmatory-test eligibility remains pending taxonomic linkage, supported phylogenies and independent-transition review; no ecological association test has been performed.

### Historical three-marker support snapshot

The initial three of the 125 planned marker trees completed inference and passed checks of source provenance, taxon membership, nonnegative finite branch lengths and resolved unrooted edge counts. Their 1,422 internal branches include 1,390 reported SH-aLRT values and 32 unreported values, retained explicitly as missing rather than zero. Support summaries and canonical unrooted split identities are archived as an incomplete diagnostic snapshot. SH-aLRT values are not bootstrap percentages or posterior probabilities. This historical snapshot has been superseded by the complete 125-marker audit described below. Model adequacy, biological interpretation of discordance and final species-tree selection remain incomplete.

### Full-proteome domain input preparation

All 526 representative proteomes were checksum-verified and streamed into a domain-search input catalogue preserving 5,815,847 protein/taxon associations. Exact sequence deduplication identified 5,713,599 unique sequences; 58,879 matched existing marker-search inputs, and 5,654,720 additional sequences were written once. These additional sequences contain 2,520,735,813 residues. This deduplication affects search work only; original protein identities, taxon associations and unresolved gene-mapping status remain available through source provenance. The additional full-proteome domain searches and merged annotation subsequently completed, as documented below. Search-partition-specific E-values must not be treated as globally calibrated values merely because gathering thresholds are shared.

## Historical marker-domain annotation and initial prediction batches

All 30,134 Pfam 38.2 profiles were searched against 58,883 unique marker sequences using HMMER gathering thresholds. Complete-result validation retained 163,650 hits of all Pfam types, including overlapping hits, and linked annotations back to 59,840 marker proteins. The 13,983 sequences with overlaps require architecture review; no-hit states are not proof of domain absence. Additional full-proteome searches subsequently completed; their results are included in the full-data procedures below.

A checksum-pinned ESMFold v1 checkpoint and frozen missing-candidate marker queue were prepared. The initial 128-sequence production chunk is restricted to complete canonical proteins no longer than 512 residues, with longer and noncanonical proteins explicitly deferred. Execution settings, source distinctions, confidence scaling and validation are documented in structure-prediction-workflow.md. At launch, actual prediction completion remained unverified; preparation and tests are not structural results.

The initial local prediction chunk subsequently completed: 128 ESMFold v1 predictions, 129 taxon-marker links, all passing independent coordinate/sequence/confidence/PAE artifact readback. The separate loading/execution audit excluded use of the uninitialized contact-regression head and reproduced saved confidence/PAE for one sequence. All 10,394 remaining eligible short-marker candidates were then resumed with identical settings and subsequently completed. These are predicted structures, not experimentally validated structures. The 64-taxon OrthoFinder core and subsequent full assignment completed; expanded family inference and its current reconciliation status are described below.

## Domain-conditioned geometry and placement uncertainty

Pfam annotations were joined to exact-sequence AFDB marker models and existing profile residue correspondences. Single-instance non-overlapping Domain-type hits with at least half HMM coverage were assessed using at least 30 pLDDT>=70 matched residues and half shared-domain coverage. Each whole-marker baseline was reproduced. Identical domain residues were compared under whole-marker and independent proper rigid fits, and local residue-distance differences were assessed under directional PAE filters. A separate analysis compared all nonadjacent cross-domain residue pairs at PAE cutoffs 5, 10 and 15 Å in both directions in both models. There were 738 domain comparisons and 375 domain-pair/taxon-pair placement assessments. These remain dependent descriptive observations; fitting improvement, high global RMSD and raw pairwise sequence differences are not branch rates, selection tests or proof of biological domain motion. Detailed filters, exclusions, sources, code and an artifact-control example are in domain-structure-comparisons.md.

## Coordinate-derived 3Di benchmark

Native Foldseek encodings (pinned commit e3fadcd07f971e864c094ac4f3a78bf4ed845e07) were extracted for the 422-model AFDB snapshot. The complete amino-acid sequence, native 3Di state string and all ten geometric features were verified; a separate scripted extraction reproduced every model. Spatial partners/features were reconstructed from original backbone coordinates, with native descriptor rounding accounted for. Invalid terminal states were masked separately because they share the ordinary native coil symbol. Confidence masks covered either the focal residue or all six residues used by each encoding, optionally requiring every directional PAE in that context to be <=10 Å.

The benchmark reproduced 858 whole-marker and 738 domain baselines and recomputed amino-acid mismatch, 3Di mismatch and geometry on identical retained sites. Three regimes yielded 4,743 accepted and 45 excluded comparison rows; the strictest retained 856 whole-marker and 717 domain comparisons. These dependent uncorrected state fractions do not establish branch lengths, physical displacement, independent evolutionary events or model suitability. Native encoding of predicted coordinates does not eliminate sequence-derived prediction circularity. Details and reproduction are in structural-alphabet-benchmark.md.

## Initial paired sequence and structural branch fitting checkpoint

Audited coordinate-derived 3Di states were projected onto the original profile alignment with identical amino-acid/structural missingness masks. All six feature residues required pLDDT ≥70 and maximum directional context PAE ≤10 Å. The 526-taxon/125-marker audit yielded 52 eligible marker alignments (4–12 taxa each, 105 taxa in their union). IQ-TREE 3.0.1 inferred LG+F+G4 marker sequence topologies, then fitted AF+G4, AF+F+G4 and LLM+G4 published structural models on each fixed topology, keeping identical sequences and invariant columns. Publisher model files were checksum-verified. All 208 fits passed the provenance and branch-correspondence audit, yielding 416 matched branches. These are conditional point estimates with frequent rare-state/near-zero warnings; no acceleration or coupling test has been established. See paired-structural-phylogenetics.md for commands, resource estimates, diagnostics and unresolved uncertainty/model checks.

## Conditional branch resampling and feature dependence

For each of the 52 covered markers, 200 paired draws at circular block lengths 1, 10 and 30 used identical sampled alignment positions for AA and 3Di, preserving taxon identities, missingness correspondence and invariant sites. Both branch lengths and gamma shape were re-estimated on the fixed original sequence topology under LG+F+G4 and the published Q.3Di.AF+G4 model. Of 31,200 draws, two were explicitly unestimable and not replaced; 62,396 fits passed full source-column, FASTA, model, topology and branch-value audit. Conditional percentile summaries cover 416 branches in each of three regimes and two alphabets. They do not include topology, model, prediction or alignment uncertainty. The feature-overlap audit found 138,415 of 474,044 shared-coordinate feature pairs separated beyond any single circular 30-column block, directly limiting interpretation of local-block intervals as calibrated uncertainty. Median interval widths were similar across block lengths, with both increases and decreases for individual branches. Reproduction and the inspected figure are in paired-structural-phylogenetics.md.

## Tree-path versus geometric benchmark

For the 52 fitted markers, leaf-to-leaf path sums were checked against tree traversal and compared with direct CA geometry on jointly observed confidence-qualified positions. Of 730 candidate taxon pairs, 722 passed ≥50-site and ≥half-of-each-taxon overlap requirements. Path intervals preserve covariance by summing branches within each joint block-30 draw before taking percentiles. Geometry includes whole-marker proper-rotation RMSD and local CA-distance changes, with a separate both-direction/both-model PAE≤10 filter. The benchmark records that tree fits use more marker information than an individual pairwise geometric subset. It does not establish physical branch displacements, geometric additivity, acceleration or independent pairwise statistics. The maximum-RMSD pair (5000823at2759, hybrid F332112 versus F4918) was diagnosed using a shared Pfam repeat region: its independent fit was 1.47 Å versus 31.12 Å whole-marker RMSD, and only 15.17% of cross-region pairs passed PAE. This is an artifact-control priority with annotation/taxon/source checks pending, not a validated innovation. See tree-path-geometry.md.


### Matched predictor–experiment controls and region-placement diagnostics

For completed exact-sequence controls, AlphaFold, ESMFold and experimental CA
coordinates were compared on identical observed residue sets at joint focal
pLDDT thresholds 0, 70 and 90. Whole-protein eligibility required at least 50
residues and half the canonical sequence. All deposited chains and models
were retained. Local pair masks required sequence separation of at least three
and distance at most 15 Å in any of the three coordinate sets, shared across
all comparisons. Superposition and distance calculations underwent independent
numerical readback. Descriptive summaries used medians across chains within
deposited models, models within entries, and entries within canonical proteins.
Paired predictor differences were calculated before hierarchical aggregation.

For the completed 12-control, 513–741-residue tier, every raw Pfam hit span was
also analyzed; overlapping annotations and repeat/family types were retained
explicitly and were not called independent domains. Regions required at least
20 observed residues and half their annotated span. Separate regional fits
were compared with regional residuals after fitting the full observed mask.
This comparison mechanically favors separate fits and is a diagnostic of
placement sensitivity, not a significance test or demonstration of motion.
ESMFold PAE was summarized within and between annotated spans over full
canonical positions, retaining both matrix directions and excluding within-hit
diagonals. Distinct overlapping hits were excluded from between-hit PAE
summaries. Native ESMFold and serialized AlphaFold focal confidence determined
those masks, which need not equal experimental-coverage masks. Quantiles were
computed in float64 and independently read back from scalar matrix values.

The illustrated Q04305 case was selected after observing large whole-protein
discrepancy and remains exploratory. Neither low local RMSD, focal confidence,
nor PAE establishes general accuracy, training independence, biological motion
or functional divergence. Experimental assembly/construct context and broader
sampling remain unresolved. Commands, complete exclusions and evidence are in
experimental-structure-workflow.md. Other prediction tiers and integrated
analyses are tracked separately and are not implied complete by this section.


### Marker-to-species edge correspondence

To expose ambiguity when marker trees cover different taxa, full-guide
bipartitions were restricted to each marker's taxon set and canonicalized as
unrooted splits. A matching marker edge was classified as a unique full edge
or a path formed by several full edges after pruning. Absent matches and
differences between the profile and MAFFT guides were retained explicitly.
No marker branch estimate was distributed across a collapsed path. Independent
repeated taxon pruning checked all projected compatibility outcomes and the
corresponding guide-edge length sums. The completed projection covers the
89-marker combined ESMFold and 124-marker expanded AlphaFold cohorts, kept
separate. The guides lack support estimates, so this analysis establishes
conditional mapping eligibility rather than a final clade assignment,
reconciliation or acceleration test. See marker-species-edge-projection.md.


## Full-data phylogenetic checkpoint, 22 September 2026

All 125 marker trees were inferred with IQ-TREE 3.0.1 on the frozen profile
matrix partitions and their documented taxon-coverage masks. The recorded
model search considered LG, WAG and JTT exchangeabilities with empirical
frequencies and gamma rate variation; each run used 1,000 SH-aLRT replicates.
A complete independent readback reconstructed 23,441,199 retained alignment
characters from the original matrix and coverage exclusions, and recovered
59,315 internal splits by graph-edge removal. Unreported support remains
missing; SH-aLRT values are not bootstrap percentages or posterior probabilities.

For each of the two homogeneous guides, each internal guide split was
restricted to each marker's actual taxon set. Splits with fewer than two taxa
on either side were uninformative. Exact marker splits and incompatible
splits (all four bipartition intersections nonempty) were classified at
descriptive SH-aLRT cutoffs of 80 and 95. All 261,500 guide/marker/edge/cutoff
rows and 2,092 edge summaries were checked. Independent exhaustive maximum
conflict searches covered five deterministic guide edges per marker/guide;
the remaining maxima were checked through their witnesses, not independently
reoptimized over every candidate. Restricted splits can represent collapsed
paths and are not independent branch replicates or gene concordance factors.
The separate role-separation diagnostic tested an exact unrooted split between
retained fungal and outgroup entries in each tree, requiring at least two
entries of each role. It does not locate a root along an edge or infer event
direction. Full commands and receipts are in
[complete marker diagnostics](complete-marker-diagnostics-20260922.md).

On the 49,027-site profile alignment, both completed guide-conditioned C20-PMSF
runs retained all 526 taxa and 1,000 bootstrap trees. Saved site-frequency
vectors, tree/report identities, replicate tip sets and branch values were
validated; empirical split frequencies were reconstructed from every replicate.
Support remains attached to explicit splits rather than transferred by node
number. SH-aLRT calculations and likelihood optimization were not independently
rerun. Guide-conditioned profiles differ, so cross-run likelihood differences
are not a direct topology preference test. The other crossed alignment/guide
runs, model adequacy, rooting and marker/taxon sensitivities remain required.
See [PMSF validation](pmsf-profile-readback-20260922.md).

## Completed local predictions and residue correspondence

Five disjoint exact-sequence ESMFold v1 cohorts contain 25,322 models: 10,522
original short proteins, 4,252 additional short proteins, 675 ecology-targeted
proteins, 5,510 proteins of length 513–768 residues and 4,363 of length
769–1,024 residues. The earlier 5,121-model partial original snapshot overlaps
the complete original cohort and is not added to these totals. Seven native
prediction configurations remain recorded separately; their differing fields
are input-receipt hash, maximum permitted sequence length and visible GPU.
This does not establish the absence of batch effects.

Original model records, exact protein sequences, prediction receipts, coordinate
hashes and directional PAE bindings were preserved in a combined inventory.
Mapping used complete-sequence identity, including reuse across identical
proteins while retaining all taxon/marker/protein associations. Independent
readback reconstructed non-gap retained Stockholm positions using cumulative
residue counts and checked matrix amino acids and confidence values against
the original prediction arrays, allowing documented PDB serialization rounding.
All 25,509 marker/protein links and 8,335,615 matrix-residue links passed. The
readback uses a shared alignment parser and does not independently infer the
alignment or validate biological coordinate accuracy.

Availability in this completed local inventory was joined with the refreshed
frozen AlphaFold catalog by exact marker, taxon, protein and sequence identity.
The union covers 56,171 of 59,840 recovered marker records (93.87%), with 722
records in both catalogs, 24,787 ESMFold-only and 30,662 AlphaFold-only records.
The 3,669 uncovered records represent 3,656 unique sequences. These availability
counts precede confidence qualification and are not whole-proteome coverage.
Absence from these catalogs is not evidence of public-database absence.

All five ESMFold cohorts completed native structural-feature extraction,
coordinate readback and directional-PAE context checks against original arrays.
Of 12,221,520 whole-protein residues, 8,306,117 (67.96%) pass native-state
validity, pLDDT >=70 at all six encoding-context residues and maximum ordered
context PAE <=10 Å. Counts precede alignment eligibility and do not measure
predictor accuracy. Model length and sampling differ across cohorts.

AA and 3Di characters were projected onto the same source alignment using
identical observation masks. Taxa required at least max(50, ceil(0.30 times
the source alignment length)) jointly observed positions; markers required
at least four eligible taxa. Only columns missing in every eligible taxon
were removed, retaining invariant columns. Independent reconstruction checked
every character and all 526-by-125 eligibility combinations. The completed
inputs contain 122 markers, 294 taxa, 44,198 retained marker columns and
6,758,598 observed paired cells. Expanded point fits and block-resampling
uncertainty remain in progress; these input counts do not update the earlier
cohort-specific evolutionary results. Refreshed AlphaFold feature/confidence
processing remains a separate active workflow. See the
[complete prediction checkpoint](completed-prediction-inventory-20260922.md)
and [expanded paired-input workflow](expanded-paired-inputs-20260922.md).

## Whole-proteome structural availability

All 5,815,847 representative proteins across 526 sampled entries were screened
against a frozen AlphaFold inventory using exact sequence hash and length.
The source-specific selection retained one model per sequence, ranked by mean
Cα pLDDT, version and model identifier. It covers 1,319,513 protein links
(22.69%) with 1,290,278 unique models across 496 entries. Identical-sequence
reuse retains each species/protein link and is not an independent observation.
The catalog excludes local ESMFold additions and precedes confidence filtering.

The producer checked selected coordinate-file hashes against existing
verification records. Independent readback reconstructed all representative
FASTA identities, model selections and per-taxon coverage from the frozen
inventory; it did not repeat coordinate-content or confidence validation.
Full source-family joins under both guide partitions were independently
reconstructed with Python sets. Model availability is not validated orthology
or evidence of statistical power. CPU Foldseek database construction is in
progress; structural clustering and remote-homology analysis are unfinished.
See [whole-proteome catalog methods](whole-proteome-structure-catalog-20260922.md).

## Full-proteome annotations and family architecture representation

The merged Pfam database retains 8,103,610 raw gathering-threshold hits across
5,713,603 unique search queries, including four marker-only queries, and all
5,815,847 representative protein links across 526 taxa. Raw scores, coordinates,
annotation types, source partitions and E-values remain available without
rescaling partition-specific significance values. Full independent checks
covered hit fields and protein identities.

Candidate within-clan competition crossed alignment versus envelope spans
with independent domain E-value versus domain bit-score ranking. Decimal
arithmetic preserved small E-values and ranking ties. Deterministic secondary
ordering did not remove primary-score uncertainty. Curated directed nesting
relationships combined with strict containment could preserve candidate nested
hits; cross-clan and missing-clan overlaps remained unresolved. Every retained
and suppressed hit, blocker and uncertainty disposition was preserved and
independently reconstructed under all four policies. This is the documented
project policy, not a claim to reproduce an unspecified PfamScan release.

Candidate architectures retain Pfam entry types, repeat multiplicity, coordinates
and alternative retained sets. Coordinate ordering and model/type signatures
are descriptive representations, not proof of biological boundaries or domain
homology. The database contains 5,767,208 alternatives summed across queries;
this is not a count of globally distinct architectures. No-hit queries remain
missing annotation evidence. Every candidate field, alternative reference,
policy uncertainty count and protein link passed independent readback. See
[candidate architecture methods](candidate-domain-architectures-20260922.md).

The full annotation database was joined through native gene IDs to both
expanded family partitions (658,183 profile-guide families and 658,522
MAFFT-guide families). Independent validation checked every native identity,
sequence/source link, family member and source crosswalk, plus database
integrity and final hashes. Family indices are guide-specific: comparisons
use full membership identities rather than matching family labels. These
partitions are not completed reconciled orthology. The profile-guide native reconciliation completed; MAFFT-guide reconciliation
is running on an isolated copy. All 525 profile HOG tables passed source
identity, species-clade containment and within-level uniqueness checks, covering
9,277,304 HOG rows and 93,034,008 assignments summed across levels. These are
not distinct-gene counts or proof of ancestral membership completeness.
Of 5,315,734 nonsingleton-family proteins, 41,129 lack root assignments;
35,187 have native misplaced-gene flags and 5,942 do not. The latter are present
in resolved trees but outside emitted root HOG parent clades. They remain
explicitly unassigned rather than being treated as biological losses. Native
classification replay across all affected families is running; event semantics
and alternative-guide validation remain outstanding. See
[HOG validation](hierarchical-orthogroup-validation-20260922.md).

Within-family architecture summaries completed across all four policies and
both full guide partitions (5,266,820 family/policy rows). They distinguish
ordered signatures, multiplicity-preserving multisets and model/type sets,
retain missing hits and ambiguity counts, and compare observed data with a
conservative subset excluding policy disagreement, rank ties, overlap, candidate
nesting and retained HMM matches below 0.70 coverage. The conservative rule is
a sensitivity choice, not proof of architecture correctness. Independent full
summary reconstruction passed. Exact-membership comparisons identified 656,507
shared families and verified identical metrics for all 2,626,028 shared
family/policy rows. Guide-specific membership remains explicit. Across either
guide, 2,023,005 proteins (34.78%) lack qualifying Pfam hits; missing annotation
cannot establish domain loss.
Domain gain/loss, fusion, rearrangement and duplication tests remain unfinished
and must incorporate supported genealogies, annotation uncertainty and
missingness. See [family-domain linkage](family-domain-bridge-20260922.md) and
[architecture variation](family-architecture-variation-20260922.md).

The completion evidence used for this update is indexed, with receipt hashes
and scoped statuses, in `metadata/methods_checkpoint_20260922_sources.json`.

The September 23 structural-coverage, confidence, architecture and HOG updates
are indexed separately in `metadata/methods_checkpoint_20260923_sources.json`.
Earlier executed analyses above retain their original cohorts and scopes.

## Completed domain structural atlas and representation census

The September 23 domain atlas retained both alignment and HMM-envelope
boundaries for the union of eligible Domain hits across the four documented
annotation policies. An interval was identified by the full source-sequence
SHA-256 and inclusive residue coordinates; hit and boundary associations
remained explicit when intervals coincided. The 594,797 model/hit pairs yielded
1,078,592 distinct intervals from 427,255 source models. Every exported atom,
residue, coordinate, occupancy and confidence value was read back against
the original CIF arrays and fragment sequence across 428 archive shards.
No interval was rejected or missing backbone atoms. The shared CIF lexical
parser and absence of biological boundary/PAE validation limit this audit.

Foldseek database construction used `--gpu 0 --threads 4
--mask-bfactor-threshold 70 --coord-store-mode 1`. Full sequence and coordinate
readback covered 168,431,396 residues. Candidate clustering used
`--alignment-type 2 --cov-mode 0 -c 0.8 -e 0.001 -s 7.5
--max-seqs 1000 --cluster-mode 0 --cluster-reassign 1`, with eight threads
and a 64 GiB split-memory limit. The exact binary hash, commands and source
hashes are pinned in `metadata/full_domain_search_database_plan.json` and
`metadata/domain_clustering_plan.json`. Database confidence masking is not
equivalent to full residue or domain confidence qualification.

The partition contained 70,537 candidate clusters, including 41,930 singletons.
Every interval occurred exactly once and each representative belonged to its
own cluster. All 594,797 alignment/envelope pairs were reconstructed from
the original boundary associations and partition: 111,002 had identical
intervals, 397,292 distinct intervals shared a cluster, and 86,503 distinct
intervals occupied different clusters. Disagreement was therefore 17.88% among
the 483,795 nonidentical pairs. Proportions were also summarized by envelope
extension length. This assesses boundary sensitivity within a joint partition;
it does not measure stability across separate clustering runs, validate native
alignment thresholds, or establish structural evolutionary change.

Domain intervals and clusters were joined to 439,217 protein records, their
taxa, and both guide-specific protein-family assignments. The 688,246 distinct
cluster/protein links count a protein once within a cluster while retaining
its original intervals and allowing membership in multiple clusters. An
independent implementation rebuilt every source identity, relation and summary
using Python sets and counters. Structural clusters were not treated as
orthogroups or proof of remote homology.

The representation census retained all 526 manifest taxa, with all six counts
(models, proteins, intervals, clusters and the two family counts) agreeing
between set-based and SQL aggregation. The selected domain catalog represented
452 fungal entries and 23 outgroups; 49 fungal entries and two outgroups had
zero representation. These zeros indicate missing evidence in the selected
AlphaFold/Pfam catalog, not biological domain absence. Alternative boundaries,
shared models, uneven coverage and shared ancestry preclude treating entries
as independent evolutionary observations. Confidence and clustering-parameter
sensitivity, direct within-domain comparisons, and phylogenetic tests of domain
and structural evolution remain unfinished.

Evidence for this completed-methods section is indexed in
`metadata/methods_domain_atlas_checkpoint_20260923_sources.json`.

The domain partition was additionally joined to the source Pfam accessions
and clans using `scripts/summarize_domain_cluster_pfam.py --plan
metadata/domain_cluster_pfam_plan.json`. Each model/hit identity and both
boundary coordinates were checked against the audited domain registry before
aggregation. The resulting 1,189,594 boundary links formed 105,519
boundary/cluster/accession/clan groups. Counts of source models, model/hit
pairs and intervals were retained separately. Source-derived Python sets
agreed with SQL aggregates, and every serialized summary row was checked.
Within the alignment-boundary view, 1,305 of 52,800 represented clusters
contained multiple Pfam accessions; within the envelope-boundary view,
1,207 of 48,658 did. Both views derive from the same joint partition and
include shared source models. These are annotation-composition counts,
not independent evidence of remote homology, functional equivalence or
evolutionary events. Candidate interpretation still requires confidence
qualification, checks for overlapping annotations and shared models, and
direct comparisons within supported genealogies.

## Direct comparisons of cross-clan domain candidates

All within-cluster Pfam accession pairs were enumerated separately for the
alignment and envelope boundary views. Source-model overlap was retained
explicitly. Models exclusive to one accession within a given pair and cluster
were used for direct comparisons; this local exclusivity does not establish
global annotation exclusivity or independent evolutionary observations.
The different-known-clan subset contained 68 boundary/cluster/accession-pair
entries representing 41 unique accession pairs, 2,669 source models and
3,559 domain intervals. Unknown-clan pairs were retained in the broader screen
but were outside this particular cross-clan comparison subset.

Candidate coordinates were extracted from the audited domain archives. Every
selected member and written PDB matched the recorded checksum. Fragment
sequence, complete C-alpha coverage, finite coordinates and mean confidence
within PDB rounding tolerance were checked. The 17,769 unique unordered
interval pairs were aligned in both input orders with the checksum-pinned
US-align 20241108 executable using `-mol prot -mm 0 -outfmt 0 -ter 2`.
All 35,538 directed comparisons completed, retaining native raw output,
alignment strings, both length-normalized TM-scores, aligned length, RMSD
and sequence identity. Both orders were retained to evaluate heuristic
alignment asymmetry. These were full-domain alignments without a PAE mask or
removal of low-confidence residues before fitting.

An independent implementation reconstructed every residue correspondence,
aligned length and sequence identity, and recomputed the least-squares RMSD
using proper rigid rotations. All RMSDs agreed with native output within
its two-decimal rounding; the maximum discrepancy was less than 0.005 Angstrom.
TM-scores were checked against native output but were not independently
reoptimized. Matched-residue confidence was calculated from rounded PDB
C-alpha pLDDT values. Input/output checksums and an independently reconstructed
pair grid established complete coverage of this candidate subset.

Descriptive summaries used the minimum TM-score across both length
normalizations and both input orders, minimum coverage and matched confidence,
maximum RMSD, and maximum score difference between input orders. Confidence
fractions of 0, 0.5, 0.8 and 0.9 were applied to both full-domain fractions and
both matched-residue fractions with pLDDT at least 70. All 68 entries remained
in each summary, including strata with no retained comparisons. A prioritization
screen of minimum TM-score at least 0.5 and minimum coverage at least 0.8
retained at least one interval comparison for 23, 23, 20 and 18 entries,
respectively. These cutoffs are descriptive choices, not significance tests
or proof of homology. Shared models, overlapping boundary views, selection
through structural clustering, annotation uncertainty and missing PAE
qualification limit biological interpretation. Functional and phylogenetic
validation remain unfinished.

Completion evidence is indexed in
`metadata/methods_cross_clan_checkpoint_20260923_sources.json`.

## Exploratory focal domain-profile sensitivity

Five full-length proteins from the PF00085/PF26973 shared-family case were
examined for annotation sensitivity before interpreting an apparent domain
transition. These proteins were selected after the structural candidate screen;
this was an exploratory case analysis. Source-sequence hashes were checked
against the archived export. PF00085.27 and PF26973.1 were extracted from the
checksum-verified Pfam 38.2 library with `hmmfetch`. HMMER 3.4 `hmmsearch` used
one worker thread and seed 42 under three conditions: standard `--cut_ga`,
`--max --cut_ga`, and diagnostic `--max -T -1000 --domT -1000`.
The first two conditions test sensitivity to acceleration filters while
retaining the curated sequence and domain thresholds. The permissive condition
records weaker matches without accepting them as annotations. Option semantics
follow the [HMMER manual](https://github.com/EddyRivasLab/hmmer/blob/master/documentation/man/hmmsearch.man.in).

Both gathering-threshold searches returned seven identical domain-hit rows.
The permissive search returned 22 rows; requiring both gathering thresholds
recovered exactly the same seven rows. Sequence/domain scores and explicit
alignment, envelope and profile coordinates were retained. E-values refer to
the five-sequence target database and were not compared directly with those
from the full-proteome searches. All 36 rows across the three conditions were
reparsed using Bio.SearchIO, including coordinate conversions and threshold
flags; output checksums were verified. This checks the recorded searches, not
biological homology or independent evolutionary events. Repeated-domain
correspondence, gene-tree support and ancestral interpretation remain pending.

Commands, complete hit records, evidence and limitations are documented in
[the focal case record](cross-clan-family-mapping-20260923.md#focal-profile-sensitivity--september-26).


## Full-cohort conditional site coupling and copy-review sensitivity

For the completed ESMFold marker cohort, we analyzed 44,198 aligned sites
from 122 markers. Each of 24 specifications regressed log(1+structural-alphabet
site rate) on log(1+amino-acid site rate), median relative solvent accessibility
(RSA, centered at 0.25), their interaction, coverage and lower-quartile
C-alpha prediction confidence, with marker intercepts. The grid crossed
three structural-alphabet models, Gamma versus FreeRate site-rate estimates,
two RSA normalizations, and inclusion versus exclusion of amino-acid entropy
and composition controls. QR least squares was independently checked against
SVD least squares and manually assembled finite-sample-corrected
marker-cluster covariance. A pseudoinverse-based precursor failed its
coefficient consistency check and was preserved as a failed run; no
tolerances were relaxed for recovery.

Inference used marker-cluster t references and BH correction across 72 focal
coefficient tests per analysis. We additionally used 2,000 whole-marker
bootstrap draws per specification and omitted every marker in turn.
Bootstrap percentile intervals are unadjusted sensitivity intervals. An
explicit paired analysis excluded marker 4986044at2759 (307 sites) because
of unresolved TFIIB/BRF1 copy assignments. Remaining covariates were identical,
and all focal omission coefficients were verified against independently
absorbed leave-one-marker-out estimates from the retained cohort. Both
analyses and all model specifications are reported; they do not constitute
independent replications. These conditional analyses do not propagate rate
or topology estimation uncertainty, fully resolve phylogenetic dependence
across markers, or remove circularity from sequence-derived predictions.
See [conditional site coupling](conditional-site-coupling.md) for provenance,
reproduction, numerical diagnostics, results and the full sensitivity figure.


## Provisional extant references for terminal duplicate pairs

For each structurally covered terminal duplicate tip pair, we inventoried
modeled nonfocal proteins in its immediate resolved-tree sister clade.
Provisional eligibility required a bifurcating parent absent from the native
duplication table and no focal-taxon genes in the sister clade. Missing models,
reported parent duplications and multifurcations were retained explicitly.
The nearest modeled reference was defined using sequence-tree path length;
all ties within 1e-12 were retained, with a lexical representative for
deterministic bookkeeping. Pair identities, eligibility and complete tied
reference sets were compared between reconciliation guides. These are
provisional extant references; unreported duplication is not evidence of
speciation, and rooting, hidden paralogy, structural quality and reference
sensitivity require additional assessment before asymmetry inference.


### Expanded AlphaFold paired-fit provenance (September 26)

Following structure recovery, paired-input comparisons identified 95 markers
with unchanged amino-acid/structural-state alignments and 30 requiring refitting.
The expanded 125-marker point-estimate collection uses the 380 audited fits
from unchanged markers and 120 audited updated fits, retaining four fit settings
per marker and all source-run identifiers. Complete table readback verified
500 selected fits and 61,893 paired branch records against the original native-run
audits. No likelihoods from different observed alphabets are directly compared.
Warnings, including near-zero branches, are retained; these point estimates
alone do not establish structural acceleration or positive selection. Separate
uncertainty analysis is required for this cohort. Reproduction and provenance
are described in [the recovery record](recovery-20260926.md#completed-expanded-alphafold-fit-collection).


## Matched prediction-source sensitivity in the expanded cohorts

We identified all taxon–marker combinations eligible in both the recovered
AlphaFold and completed ESMFold paired-input datasets. Each independently
filtered alignment was mapped back to the original marker-alignment columns;
comparisons required observed, matching amino acids under both source-specific
confidence masks. The shared set contained 673 taxon–marker combinations across
21 taxa and 78 markers. Complete native-encoding sequences were then compared
by exact strings, lengths and sequence hashes. All shared combinations had
identical complete sequences, resolving to 643 distinct source model pairs.
Model versions, coordinate and encoding checksums, and available prediction
provenance were retained. Multiple taxon–marker links to one model pair were
collapsed for full-model sensitivity summaries, with the source links preserved.

For each distinct pair, we compared valid native structural-alphabet features
across the full encoded proteins. Joint masks used minimum feature-context
pLDDT cutoffs of 0, 70 and 90 and maximum directional feature-context PAE cutoffs
of unfiltered, 5, 10 and 15. This produced 12 alternatives for every pair.
A pair qualified for a threshold summary if at least 50 residues and half of
the full sequence were jointly retained. All failed coverage dispositions were
kept. We reported pooled residue disagreement, pair-level summaries and the
retained pair/residue denominators. Structural-state disagreement was also
partitioned by whether the native partner-residue index agreed between models.
This partition describes association with partner context, not a causal
attribution of disagreement.

To distinguish changes in protein membership from within-protein filtering,
we repeated all threshold summaries on the intersection of pairs qualifying
under every alternative (119 pairs). Protein membership was thus fixed, but
residue membership remained threshold-dependent. We retained all pair-specific
threshold differences; no significance test or calibrated error estimate was
constructed from these descriptive controls. Independent reconstruction checked
all source-overlap joins, complete sequences, full-model confidence selections,
partner/state counts, summaries and fixed-cohort membership.

Local model provenance was concentrated in the ecology acquisition batch
(635 of 643 pairs), with five original-marker and three longer-marker models.
The three saved configurations shared the reviewed checkpoint/software,
precision, recycle, seed, chunk-size and inference-mode settings. This does not
make the overlap representative or establish complete configuration equivalence
with AlphaFold. Confidence selection, nonlocal structural features, shared
proteins and acquisition ascertainment limit inference. These comparisons assess
source sensitivity on matched sequences; they do not measure experimental
accuracy or evolutionary substitutions. Reproduction, receipts and the
confidence/coverage figure are in [prediction-source controls](prediction-source-controls.md#qualified-full-cohort-overlap-diagnostic-september-26).

## Matched domain contrasts and phylogenetic working models (September 27)

The matched analysis retains 192 combinations of domain boundaries, confidence
masks, cohort, coverage screen and alignment input order, each with 432
species-guide, annotation-policy and matching-scenario strata. The resulting
82,944 settings are represented by an audited map to 28,808 exact record inputs.
Identity requires agreement of both the complete numeric matrix and the ordered
target, background, family-component and species-pattern identities; numeric
similarity alone is insufficient for reuse. Each input is fitted under five
species-tree alternatives, yielding 144,040 unique fits and 414,720 logical
setting/tree results. These alternatives are correlated sensitivity analyses,
not independent replicates. The full fitting and validation workflow remains
in progress.

The response is target-minus-background Cα RMSD, averaged equally over eligible
domains within each selected record. The four adjustment variables are the
corresponding exact sequence-identity difference, original-interval coverage
difference, log aligned-length ratio and joint high-confidence residue-fraction
difference. The fixed design contains an intercept and varying covariates scaled
by their population standard deviations without centering. Zero constant
columns are explicitly omitted; a nonzero constant column stops that fit for
review. Thus the intercept refers to zero covariate differences, rather than
the sample-average covariate combination. Full design checks found no remaining
rank deficiency after constant-column removal. The confidence-difference column
was constant in 41,472 settings. These checks do not establish adequate overlap
or remove sequence/prediction confounding.

For record i, the working species effect is the focal-species effect minus one
half of each background-endpoint species effect. We combine weights when taxa
coincide. This additive endpoint assumption defines a zero-sum contrast matrix
W; it does not assume that structural distances themselves are additive. For
each tree, D contains patristic substitution distances, and the centered kernel
is C = −HDH/2, with H the centering matrix. D is not squared. The species term
is W C Wᵀ. Because the rows of W sum to zero, it is invariant to root placement
for the same unrooted branch-length distances. It is neither a dated covariance
nor an estimate of structural change per unit time.

The working response covariance is

\[
\operatorname{Var}(y)=s^2\{I+r_b ZZ^T+r_f GG^T+r_p FF^T\},
\]

where Z indicates shared backgrounds, G indicates family components joined by
shared model identities, and F Fᵀ equals the selected species-contrast
covariance. Backgrounds must nest within family components. The three ratios
are nonnegative and the overall residual multiplier s² is positive. The full
species-pattern design has 4,568 rows and rank 242; compact factors retain this
estimable space without adding diagonal jitter. Independent tree-edge features
verified all 104,333,120 pattern-pair/tree covariance entries. Numerical agreement
establishes implementation of this working covariance, not its biological
adequacy.

We profile s² by restricted likelihood and optimize log(1 + ratio), bounded by
a variance ratio of 10,000. The optimizer retains all eight zero-component
faces: the exact all-zero face and three starts (ratios 0.05, 1 and 20) for
each other face, totaling 22 attempts. Saved candidates include failures,
objectives, convergence messages and boundary flags. A separate direct evaluator
checks candidate likelihoods and conditional coefficient calculations. An
analytic-gradient follow-up assessed stationarity where coarse finite-difference
diagnostics were inconclusive. All 658 targeted fits passed the recorded
numerical criteria after refinement and full output readback. Selected estimates
were integrated into 144,040 unique-fit rows and 414,720 setting rows; original
estimates and flags remain preserved. The completed integration and its scope
are recorded in `metadata/matched_refinement_integration_completed_20260928_v2.json`.
Passing numerical checks does not prove a global optimum; conditional covariance
matrices are not final calibrated uncertainty estimates.

A separate complete-input diagnostic tests whether the equal-covariate reference
is within the observed joint convex hull. Each covariate is divided by its
maximum absolute value (one for an all-zero column), without centering. A linear
program minimizes the largest absolute coordinate of a convex combination of
records. Saved nonnegative weights and separating vectors are checked against
the original matrices. Numerical support uses a maximum absolute barycenter
coordinate of 1e-8; separation requires a minimum projection greater than 1e-7
with a direction of L1 norm at most one within tolerance. Solver failures,
invalid certificates and unresolved boundary cases remain explicit. Marginal
range overlap alone does not pass this joint check, and convex-hull support does
not establish dense local sampling or causal comparability.

Production details and source-bound validation records are linked in
[full matched models](full-matched-working-models-20260927.md),
[species covariance factors](matched-species-covariance-factors-20260927.md),
[design diagnostics](matched-domain-design-diagnostics-20260927.md) and
[joint support](joint-covariate-support-20260927.md). Covariance adequacy,
nonlinear sequence adjustment, calibrated uncertainty and multiplicity treatment
remain required before interpreting a duplication-associated effect.

## Local coordinate checks at shared functional positions (September 27)

For all 150 exact protein positions jointly observed in the qualified AlphaFold
and ESMFold functional-annotation datasets, we compared four corresponding
Cα subsets: the focal residue and its sequence neighbors, the AlphaFold-selected
partner neighborhood, the ESMFold-selected partner neighborhood, and their
union. Each context includes the focal neighborhood; repeated positions in a
context receive only one equal weight. Correspondence uses exact full protein
sequences and native residue coordinates, not a new structural alignment.
All 600 comparisons were retained. Raw coordinate hashes and residue sequences
were checked for 186 source models.

Proper-rotation least-squares RMSD was checked against a quaternion eigensystem;
intracontext pair-distance differences were independently calculated by scalar
and vector operations. We retained confidence values and numerical geometry
status without additional confidence filtering. Summaries distinguish agreement
in the structural-alphabet state from agreement in the descriptor's selected
partner residue. The four groups contain 118, nine, 11 and 12 positions for
same-state/same-partner, same-state/different-partner,
different-state/same-partner and different-state/different-partner, respectively.
Each position recurs in all four contexts; these measurements are dependent.

The [four-panel figure](figures/functional_predictor_geometry_20260927.pdf)
shows all observations and descriptive medians. These matched-sequence controls
assess prediction-source sensitivity in local backbone geometry. They do not
measure experimental accuracy, side-chain or pocket agreement, whole-protein
displacement, ancestral changes or independent evolutionary events. Annotation
correspondences remain hypotheses about functional positions. Source records,
reproduction scripts and limitations are described in the
[functional-site workflow](functional-site-workflow.md).

## Nonlinear identity adjustment and ordinary-ML comparisons (September 27)

To assess sensitivity to linear sequence adjustment, we constructed quadratic
and cubic identity contrasts on the same eligible domains as the original
matched measurements. For matched record i, eligible domains D_i and power
p ∈ {1,2,3}, the covariate is

\[
 c_{ip}=|D_i|^{-1}\sum_{d\in D_i}
 \left(I_{id,\mathrm{target}}^p-I_{id,\mathrm{background}}^p\right).
\]

The transform precedes averaging. Squaring the mean identity difference, or
subtracting squared mean identities, would define different covariates. The
response remains target-minus-background structural RMSD; original coverage,
log aligned-length ratio and confidence-fraction contrasts remain in the design.
Each higher degree adds one identity term to the preceding design. Constant-zero
columns are removed explicitly; a nonzero constant stops the fit for review
because it changes the interpretation of the zero-covariate intercept. Active
covariates are divided by their population standard deviations without centering.
Saved coefficients and their conditional covariance are transformed back to
original covariate units.

Separate parsers verified all 9,983,040 configuration rows across 192 settings.
The expanded designs passed full independent rank and conditioning readback for
165,888 setting/degree combinations. Exact numeric bytes and ordered observation
identities define reusable inputs. All 82,944 setting labels link to 28,808
unique linear/quadratic/cubic triplets, preserving response, observations,
background reuse, family components and species patterns. Sensitivity settings
are repeated analyses, not independent biological observations.

For fixed-effect comparisons we fit all three designs under ordinary Gaussian
maximum likelihood with the same working covariance R described above. At fixed
variance ratios, generalized least squares gives residual quadratic q; the
profiled variance is q/n and the negative log likelihood is

\[
 -\ell=\tfrac12\{\log|R|+n[1+\log(2\pi q/n)]\}.
\]

This objective uses n, rather than n minus the number of fixed coefficients,
and excludes the restricted-likelihood fixed-design determinant. The earlier
REML fits remain separate and are not used as linear-reference likelihoods.
Ordinary-ML optimization uses analytic scores, the same eight component-boundary
faces and 22 total attempts, retaining failed attempts and bound contacts.
Every candidate objective is checked through a direct residual calculation.
Synthetic dense checks validate implementation but do not establish global
optimization or the biological adequacy of this covariance model.

The full ordinary-ML run covers 86,424 inputs across five tree choices, totaling
432,120 fit dispositions. Production and full output audit are in progress.
A downstream export will report unrounded differences in fitted log likelihood
for linear-to-quadratic, quadratic-to-cubic and linear-to-cubic comparisons.
A worse more-flexible fit beyond 1e-7 plus 1e-9 times the larger absolute objective
is flagged for numerical review; all negative gains remain recorded. Failed
fits and optimization-review flags remain visible in every setting/tree row.
No chi-square calibration, p-value or biological model preference follows from
this export. Joint support of the expanded zero-covariate reference is assessed
separately, retaining unresolved cases and any separately verified nonnegative
weight certificates. Predictive adequacy, uncertainty and multiplicity treatment
remain required before interpretation.

Implementation and source bindings are documented in
[nonlinear inputs](nonlinear-identity-contrasts-20260927.md) and
[ordinary likelihood](matched-ordinary-likelihood-20260927.md).

### Numerical refinement and local gradient checks (September 29)

A frozen prefix of 75,205 completed ordinary-ML fits contained 601 numerical
review flags. This prefix is an interim execution snapshot, not a random sample
or a final failure-rate estimate. All 601 source fits, including their original
failed attempts, remain preserved. Their 449 distinct numerical inputs were
reconstructed and independently checked against the original ordered records,
exact matrix signatures and covariance indices before separate refinement.
The original full model run continues independently.

Refinement retains 24 candidate records: the original unoptimized point, the
22 original face/start combinations, and an additional full-face optimization
initialized at the original parameters. Variance ratios are parameterized as
log(1+ratio), with the original maximum ratio of 10,000 retained. L-BFGS-B uses
analytic gradients, at most 1,000 iterations, objective tolerance 1e-14,
gradient tolerance 1e-8 and at most 80 line-search steps. Selection minimizes
the recorded objective; successful termination breaks only exact objective ties.
A lower objective from a failed attempt is never silently discarded.

Each candidate objective is recalculated with the direct residual likelihood.
The selected candidate must have successful optimizer termination, projected
infinity-norm gradient at most 1e-3, no upper-bound contact, full-face starting
objectives agreeing within 1e-5, and an objective no worse than the original
point. Every check is saved separately. At a lower bound the projected gradient
retains only negative derivatives; at an upper bound it retains only positive
derivatives, using the implementation's 1e-7 boundary neighborhood. Upper-bound
contact within 1e-6 remains a separate flag. A queued output audit replays every
candidate and classification; it shares the validated numerical libraries and
is not an independent mathematical implementation.

Two interim cases retained only a gradient flag despite relative-objective
convergence. Direct likelihood differences confirmed their gradients above
1e-3. Local correction diagnostics use symmetrized Hessians obtained by
differencing analytic gradients at steps 1e-5 and 1e-6. Both Hessians must be
positive definite on the adjusted coordinates. Each Newton proposal must
preserve or improve the direct objective and pass analytic and direct gradient
checks. Direct checks use steps 1e-6 and 1e-7. Interior proposals use centered
differences for all parameters. Lower-face proposals keep exact-zero variance
components fixed, use centered differences for free parameters, and second-order
forward differences for fixed components. Their analytic and direct one-sided
derivatives must be nonnegative, while free-gradient norms must be at most
1e-3. Both Hessian-step proposals must pass; checks are not chosen by which
step gives the preferred answer. Upper-bound cases and unsupported faces remain
explicitly inapplicable.

These are local numerical proposals. Original optimizer failures, disagreement
between starts and other review flags cannot be resolved by this procedure.
Single-case direct checks and synthetic fixtures have passed. Complete snapshot
refinement processed all 601 fits, with 598 passing and three gradient-only
flags retained. Full output replay subsequently passed all 14,424 candidate likelihoods. Both
proposal stages and their complete readbacks also passed, yielding two checked
local proposals for each of the three retained gradient cases. A separate selected
snapshot now passes readback for all 601 choices, including coefficients, residual
scale and conditional coefficient covariance recomputed by a separate normal-equation
implementation. Maximum likelihood discrepancy was 1.82e-12; the covariance solver
and analytic-gradient library are shared. Production integration remains pending.
No original production fit has been replaced. Completion
and validation of the full model grid remain required before using corrected
values in model comparisons. Numerical stationarity does not establish variance
component identifiability, a global optimum, interval calibration, prediction
source independence or a biological effect.

Implementations, pinned inputs, separate checks and current stage status are
linked in [numerical refinement records](nonlinear-identity-contrasts-20260927.md).

## Conditional ancestral amino-acid inference in exploratory families

We analyzed 13 exploratory families using 26 whole-protein alignments (MAFFT
and FAMSA for each family). The input union contained 1,025 proteins from 417
taxon entries and 482,007 ungapped residues; the largest family contained 622
proteins. All proteins and alignment columns were retained. These families
are mechanistic candidates, not a random sample of the fungal proteome, and
selection of a family does not establish an evolutionary effect.

For each alignment, IQ-TREE 3.0.1 fitted LG, WAG and JTT exchangeabilities with
empirical amino-acid frequencies and four gamma rate categories, holding the
unrooted topology fixed while optimizing branch lengths and gamma shape.
The initial design comprised 78 fits. We assessed optimization sensitivity
using three paired gamma-shape/branch-scale starts, (0.05, 0.5), (0.5, 1) and
(2, 2), at each of two lower bounds on gamma shape (0.02 and 0.005). These
468 alternate fits used an optimization epsilon of 1e-6. For each of the 156
alignment/model/bound groups, the best eligible audited baseline or alternate
fit initialized an additional refinement at epsilon 1e-8. All warnings,
short branches and alternative optimizer outcomes were retained.

Independent readback checked every fitted report and checkpoint, tip and edge
sets, empirical frequencies, gamma-category rates and likelihood. Frequency
reconstruction followed the version-specific eight-iteration allocation of
unknown observations. Likelihoods were recomputed by scaled pruning, with an
absolute agreement threshold of 0.001. Passing this check establishes numerical
agreement for a saved fit; it does not establish a global optimum or model
adequacy.

We computed all 20 conditional amino-acid probabilities at three identifiable
candidate internal vertices, using original labeled trees to derive partitions
of retained tips around each vertex and matching those partitions to the fitted
tree. The original degree-two root was suppressed by the unrooted model and
was recorded as having an unidentifiable position, rather than assigned to a
nearby fitted vertex. Gaps and X were treated as unknown residue observations;
this calculation did not infer deletions or an ancestral ungapped sequence.
Source alignments remained unchanged.

The probability implementation used rerooted pruning, repeated-column
compression, and spectral transition matrices, with direct matrix
exponentials for very short transitions. A separate fixed-root inside/outside
implementation used full columns and direct matrix exponentials, independently
constructing the rate matrix from the same empirical exchangeabilities and
saved parameters. An analytically enumerated two-internal-node example checked
the independent implementation. Full readback required absolute agreement
within 1e-8 for probabilities and 1e-7 for site log likelihoods. All 8,708,760
probabilities in the 156 refined fits passed this readback.

Same-coordinate sensitivity summaries compare the most probable amino acid,
its support, and total variation between the complete 20-state probability
distributions. Opposing most-probable states with at least 0.9 probability in
both fits are counted explicitly. Comparisons require matching alignment
hashes and candidate-vertex partitions. Model, bound and node comparisons are
dependent diagnostics and are not counted as independent evolutionary changes
or used as posterior ensemble weights. The refined/baseline and bound
comparisons and all alternate-start comparisons are complete. The independent
audit checked 26,126,280 alternate-start probability values. Across 1,306,314
alternate/baseline node-site comparisons, 128 most-probable states changed;
none had opposing calls supported at least 0.9 in both fits. The largest total
variation was 0.298073. A separate descriptive diagnostic joined the ten
retained per-fit maxima to observed alignment characters, descendant sets,
fit likelihoods and refined probabilities. These are not the ten largest
individual sites globally and do not establish a cause for sensitivity.

Joint insertion/deletion uncertainty is unresolved. Independent gap-character
marginals can violate mutual-exclusion constraints, and imposing compatibility
after fitting does not supply a fitted joint evolutionary model. Historian
checks revealed sensitivity of the largest-family candidate outputs to minimum
branch length. BAli-Phy initialization and short-chain diagnostics assess model
execution and sample identity, not posterior convergence. Final ancestral
sequences, ancestral structure ensembles and functional conclusions have not
been qualified by these conditional calculations.

Inputs, commands, numeric checks and completion records are linked in
[ancestral case evidence](ancestral-case-inputs.md), with separate
[Historian](historian-method-assessment-20260927.md) and
[BAli-Phy](baliphy-method-assessment-20260927.md) assessments. The
[refinement sensitivity figure](figures/whole_refinement_sensitivity_20260927.pdf)
is generated from the complete audited comparison table.


Refined whole-protein and domain ancestral marginals were compared at exact
shared projected residue-coordinate signatures across all alignment-method
pairs and both domain boundary definitions. We matched original candidate
nodes and required whole-tree descendant sets, restricted to domain-bearing
proteins, to equal domain-tree candidate sets. The three substitution models
and two parameter bounds yielded 624 fit pairs and 1,872 non-root comparisons.
Unmatched alignment columns were retained as explicit counts and in the full
source coordinate union. Most-probable residue agreement, opposing calls
supported at least 0.9 in both fits, and total variation were computed on
matched columns and independently read back against saved arrays. These
comparisons confound sequence context, membership and parameter estimation;
they are descriptive sensitivity analyses, not tests isolating context effects.

## Joint alignment and ancestral sampling: running analysis

We initiated BAli-Phy 4.3 sampling for the same 13 exploratory families to
address alignment and insertion/deletion uncertainty that is absent from
fixed-alignment amino-acid marginals. Runtime input checks reduced 324 labeled
configurations to 135 distinct combinations of ordered ungapped sequences and
fixed input trees. Every original configuration remains linked to its effective
input. Labels with identical effective inputs share computations and are not
treated as independent chains or evidence for agreement between aligners.

The substitution model uses LG exchangeabilities, amino-acid frequencies with
the installed symmetric Dirichlet prior of concentration one per amino acid,
and four gamma rate categories. Three explicit priors on gamma shape are
analyzed separately: `LogLaplace(6,2)`, `LogLaplace(0,1)` and
`LogLaplace(0,2)`. Under the audited installed implementation these have
medians approximately 403.43, 1 and 1, respectively; the first is the installed
default. The RS07 insertion/deletion model uses the installed priors
`LogLaplace(-4,0.707)` for rate and `ShiftedExponential(10,1)` for mean length.
The generated model programs make these priors explicit. Frequencies, gamma
shape and indel parameters are sampled; tree topology, branch lengths and
the scale of one remain fixed. Indel rates use the constant-rate configuration.
Consequently these chains do not propagate root, topology or branch-length
uncertainty, and their results must retain the input-tree sensitivity labels.
Prior settings are not pooled using arbitrary ensemble weights.

All 405 input/prior initializations passed generated-model and initial-score
checks. A separate fixed-parameter, 20-iteration grid passed saved-alignment
and node-identity checks for all 324 original configurations. Those short
runs supplied capacity and integrity evidence, not posterior samples for
biological interpretation. The subsequent full sampling grid comprises four
unique seeds for each of the 405 input/prior combinations, or 1,620 chains.
Its first horizon is 1,000 iterations; reaching that horizon is not a
convergence criterion. Chain identities, seeds, model sources, input hashes,
commands and resource estimates are stored in versioned manifests.

Each chain exports its labeled runtime tree from the same process immediately
before creating its MCMC state. On successful exit, readback compares clades
and branch lengths with the input tree, checks that observed tip sequences
are preserved, and maps candidate ancestors by descendant sets. Every saved
alignment is checked; process-local internal names are not matched directly
between independent chains. Scalar logs include iterations 0 through 1,000;
ancestral alignments are saved every ten iterations. Interrupted attempts
are preserved, and recovery starts a fresh attempt rather than appending
samples to a previous chain. This process-recovery path was tested with the
installed binary; native MCMC checkpoint continuation was not established.

Queued diagnostics evaluate complete four-chain groups at both 25% and 50%
burn-in cutoffs, retaining both results. Scalar screening uses ArviZ 0.22.0
rank-normalized split/folded R-hat, bulk and tail effective sample sizes, and
mean Monte Carlo error. Initial screening thresholds are R-hat below 1.01 and
bulk and tail effective sample sizes of at least 400 across the four chains.
Nonfinite, insufficient and constant stochastic traces receive explicit review
flags. Separate candidate-length screens preserve the saved iteration schedule
and stable source-node identities; they retain 75 or 50 observations per chain
at these cutoffs and do not inherit scalar-log effective sample sizes.

The sampling and diagnostic stages remain in progress. Scalar and length
screens alone will not qualify ancestral sequences: positional amino-acid
states, homology/alignment mixing, prior/root sensitivity and model adequacy
remain additional requirements. Failed or capped chains remain in the
denominator, with further sampling or methodological review required where
checks fail. No ancestral structural ensemble or functional conclusion has
yet been accepted from this joint-sampling analysis.

Exact implementation and evidence are linked in the
[BAli-Phy assessment](baliphy-method-assessment-20260927.md). The
[run plan](../metadata/baliphy_independent_chain_plan_20260927.json),
[scalar diagnostic plan](../metadata/baliphy_independent_chain_diagnostic_plan_20260927.json)
and [length diagnostic plan](../metadata/baliphy_candidate_length_diagnostic_plan_20260927.json)
define the current execution, while the
[prior audit](../metadata/baliphy_prior_definition_audit_20260927.json) binds
the distribution definitions to installed source files.

## Whole-protein sequence–structure contrasts (September 29)

This analysis extends the matched domain workflow to whole-protein structural
comparisons. Input construction, descriptive summaries and complete model-design verification
are complete; the exact-input inventory is in progress. Whole-protein
effect fitting, uncertainty calibration and biological interpretation have not
been completed. These inputs use the original matching snapshot and must not
be presented as comparisons across every model in the expanded September 28
atlas.

We retained the original 2,786,912 target/control selections across 432
guide, annotation-policy and matching-scenario strata. Their 52,675 distinct
target/control pairs were linked to canonical model/version identities and
audited structural measurements. We carried full-protein and pLDDT≥70 masks
and both alignment input orders for each member of the contrast, yielding
421,400 pair/mask/order rows. Input-order alternatives and repeated selections
are sensitivity settings, not independent biological replicates. Identical-model
cases remained explicitly unmeasured rather than receiving an assumed zero
structural distance. Missing measurements and degenerate superpositions
remained explicit exclusions.

Structural eligibility was evaluated at minimum aligned lengths of 30 or 50
residues and minimum coverage of 0.5, 0.7 or 0.9 in both original proteins.
Coverage denominators were the full, original protein lengths, including for
confidence-masked comparisons. The screen required usable measurements in both
alignment orders. We retained both mask-specific cohorts and the intersection
passing both masks, without rematching after filtering. The resulting grid
comprised 96 settings: two masks, two cohort definitions, six screens, and four
target/background order combinations. All screens were retained rather than
choosing one from observed structural outcomes. Attrition and post-screen
balance were evaluated relative to the original matching records; stricter
coverage did not automatically imply better covariate balance.

For each retained record, the RMSD response is target RMSD minus background
RMSD. The second response is the difference in divergence defined as one minus
the mean of the two endpoint-normalized native TM scores. Equivalently, it is
the background mean TM score minus the target mean TM score. These responses
compare separately fitted structural cores; they do not guarantee identical
homologous residue subsets. TM-score divergence is a structural-similarity
summary, not physical displacement. Whole-protein changes can also reflect
domain arrangement and therefore require comparison with the domain analyses.

Sequence adjustment uses the separately retained target and background
gene-tree distances, t and b. The five working fixed-effect designs contain
linear, quadratic or cubic differences of powers (t^k − b^k for k=1,…,d),
the log-distance difference log(t) − log(b), or the aligned exact
sequence-identity difference. The polynomial designs do not use powers of
(t−b). Log-distance designs require both distances to be positive; zero-distance
records are excluded with explicit denominators and no pseudocount. Gene-tree
distances are relative sequence divergence, not dated rates. Each design also
includes target-minus-background minimum original-protein coverage, the log
aligned-length ratio, and the difference in jointly high-confidence residue
fractions.

Before fitting, each of the 207,360 design/stratum combinations is assessed for
constant columns, nonzero constants, numerical rank, residual degrees of freedom,
conditioning, marginal support for zero differences, and reuse of backgrounds.
Rank diagnostics use centered and population-standardized varying columns,
with SVD cross-checked against pivoted QR. This centering belongs to rank
diagnostics; it does not redefine the zero-difference reference for subsequent
fitting. Marginal support is not joint convex-hull support. The completed
producer reported full rank and positive residual degrees of freedom in every
design; independent numerical reconstruction subsequently verified all 207,360
designs across all 96 settings.

Descriptive summaries give record-weighted means, means with equal weight per
family, and means with equal weight per focal taxon. Positive-distance record,
family and taxon denominators are recorded separately for log contrasts. These
three weighting schemes are descriptive summaries and must not be conflated
with a fitted mixed-model weighting rule. Independent reconstruction passed
all 96 settings and 1,492,992 weighted-mean cells. A historical audit mismatch
did not recur in the full replay; a separate compensated-sum check passed all
15,552 cells in the affected setting. The original failure and unresolved cause
remain documented; neither input data nor acceptance tolerances were changed.

For planned working models, all 52,675 pairs are bound to the shared-entity
family components and five species-contrast covariance factors defined above.
The focal-minus-background-endpoint contrast uses the same W construction,
including combined weights when taxa coincide. These are nuisance covariance
assumptions, not estimates of a structural evolutionary process. The input
inventory retains all 414,720 design/outcome settings and permits reuse only
when ordered pair identities and numerical input bytes agree exactly, retaining
column/outcome and covariance definitions. Signed zero is normalized; no
approximate merging or record subsampling is used. Five tree alternatives
remain separate. Completed fitting, covariance adequacy, joint support,
prediction-source controls, uncertainty calibration and multiplicity treatment
are still needed before interpreting excess structural change conditional on
sequence divergence.

Scripts, source bindings, numerical definitions and validation receipts are
linked in [whole-protein matching and validation](duplication-control-balance-20260927.md).
The [design-producer completion record](../metadata/whole_protein_model_designs_v2_producer_completed_20260929.json)
and [input-inventory plan](../metadata/whole_protein_model_input_inventory_plan_20260928.json)
distinguish completed preparation from pending model inference.


The expanded September 28 duplication queue contains 134,812 distinct model/version
pairs involving 276,682 models. Exported C-alpha sequences, coordinates and
confidence values for all 100,073,779 residues were reconstructed independently
from the frozen source mmCIF atom rows. The audit shares the lexical CIF parser
with the producer but uses separate extraction logic. All 277 shard proofs and
aggregate counts passed; there were no rejected model dispositions. Source and
proof bindings are recorded in the
[completed audit](../metadata/duplication_coordinate_readback_completed_20260929.json).
Full and pLDDT70 C-alpha alignment inputs retain original residue positions;
their full source/PDB readback and expanded alignment/input-order checks have
completed, as described below. Coordinate auditing provides no PAE
qualification or biological effect estimate.

## Expanded structural input-order checks (September 30)

The expanded cohort completed 539,248 pair/order/mask dispositions, retaining
both full-protein and pLDDT70 inputs and both orders for every distinct pair.
All 501,324 successful alignments underwent independent coordinate geometry
reconstruction. Numerical usability required at least three matched residues,
the unchanged RMSD agreement tolerance, and a unique optimal proper rotation
under the audited numerical criterion. Excluded directions retained their
native status and overlapping reasons; their numerical metrics remained blank.
This is a numerical screen, not confidence or coverage qualification.

After endpoint-normalizing both orders, the complete order summary was
independently reconstructed. Native aligned-residue correspondences for every
usable two-order pair/mask comparison were then replayed with separate
cumulative-index arrays and compared as sets. We recorded identical mappings,
changed mappings of equal size, changed counts, set Jaccard similarity and
absolute order differences in all eight metrics. Independent sorted linear
interpolation reproduced every metric quantile. Full comparisons and the
same-pair common-mask cohort remained separate; no order was favored or
averaged. RMSD order differences can involve different residue correspondences
and are not uncertainty bounds or evolutionary displacement. Full evidence and
descriptive figures are linked in
[expanded input-order sensitivity](duplication-sampling-coverage-20260926.md#expanded-input-order-sensitivity-completed-september-30).

## Marker topology and supported-conflict sensitivity (September 30)

All 125 MAFFT marker inputs were reconstructed from the frozen source matrix
and unchanged coverage rule. Independent graph-edge removal recovered every
internal support split. Profile/MAFFT marker-tree comparisons restricted split
sets to shared taxa, removed trivial projections and collapsed duplicate
projected splits, while retaining zero-length internal edges. Independent
DendroPy pruning checked the complete comparison. Restriction to shared tips
is not refitting: inference effects of original taxon membership remain.

For every full homogeneous guide-tree edge, each marker split was restricted
to its retained taxa. We retained exact matching SH-aLRT support and the
highest support among incompatible marker splits, requiring all four
bipartition intersections to be nonempty for conflict. Cutoffs 80 and 95
produced explicit concordant, conflicting, unresolved and uninformative-coverage
statuses; missing support was never interpreted as zero. Independent Python
set searches checked every maximum separately from the producer's bit-vector
implementation, and verified that both cutoff rows shared the same maximum.
All exact supports, witnesses, complete grids and aggregate counts were checked.

An independent outer one-to-one merge reconstructed all cross-method
classifications and transitions from the audited source tables. Every guide-edge
sensitivity summary was retained, together with separate strata for the 95
markers with identical original tips and the 30 with different memberships.
SH-aLRT cutoffs are descriptive settings, not bootstrap or posterior
probabilities. Rows share markers, taxa, guide edges and possibly collapsed
paths; these counts are neither independent replicates nor gene concordance
factors. No alignment method was selected as correct, and no biological cause
or accepted species topology follows from disagreement. Species-tree model,
root and calibration uncertainty remains required for downstream structural
branch tests. See the [completed full workflow](marker-alignment-topology-comparison-20260927.md).

## Expanded coverage and confidence-mask attrition (September 30)

We applied the existing six combinations of 30/50 aligned residues and
50/70/90% coverage to every expanded pair/mask comparison. Both input orders
must be numerically usable and meet the residue and coverage thresholds.
Coverage denominators are the original full-protein lengths from the frozen,
sequence-exact model queue, including for pLDDT70 inputs. The producer uses
exact integer rational comparisons; separate readback uses decimal ceiling
cutoffs and SQL model/length joins. All order-specific numerical exclusions,
screen exclusions and both-mask intersections remain attached. Confidence
masking can alter correspondence sets, so masked and full passing cohorts
are not assumed to be nested.

All decisions were projected onto the complete previously audited terminal
singleton-side candidate event ledger, retaining both guides and every source
field. Model/version associations follow gene identity even when the source
ledger and reviewed queue reverse gene positions. Same-model comparisons stay
unmeasured; missing and unresolved records retain explicit exclusions. No
event is replaced by zero structural change. Independent SQL event joins and
recomputed pair decisions verified every event screen and aggregate. Every
guide/taxon/mask/screen cell for all 526 working entries was retained, including
zero-pass cells. The full terminal candidate ledger and full manifest remain
the denominators in descriptive event and taxon summaries; guides and repeated
events are not treated as independent replicates. This screen does not qualify
the species sampling as unique taxa, establish confidence calibration or
domain orientation, correct ascertainment, or test a biological duplication
effect. Full evidence and figures are in
[expanded coverage screening](duplication-sampling-coverage-20260926.md#expanded-original-length-screening-completed-september-30).


### Expanded duplication/control integration (September 30)

The September 28 primary duplicate queue was annotated in full against the
independently verified expanded Pfam registry: 276,682 model/version identities,
134,812 distinct model pairs and four alternatives per pair. Independent segment
and policy-membership joins reconstructed every annotation and all 539,248 pair
categories/candidate single-copy domain-match records. Repeated accessions remain
explicit; a retained ineligible repeat cannot make an eligible occurrence unique.
Identical-model event endpoints remain in the model inventory. Missing annotation
is unknown, rather than evidence of a gain or loss. These controls do not establish
structural domain boundaries, homology, evolutionary architecture events or stable
interdomain orientation. Expanded sibling-reference comparisons remain separate.
Completion evidence is `metadata/expanded_duplication_domain_controls_completed_20260930.json`.

The background refresh retains every bifurcating terminal sister pair in the
unchanged two complete resolved-tree sets and rejoins the frozen expanded structure
bridge. Independent native-tree, guide-union, native reciprocal orthology,
model/metadata and raw-registry checks are running or queued. Computational
qualification preserves all three prior native-orthology/parent sensitivities,
all four annotation policies and all sequence-distance/focal-taxon ranges over
the full 283,409 modeled duplicate target links. Future artifacts are checked
against complete verified source manifests before and after architecture counts,
not accepted merely because a dependent process exited. Full expanded candidate
graphs, all 54 response-independent matching scenarios and coordinate measurements
remain outstanding. The existing completed selections and ongoing fitted models
retain their older frozen source scope. No rematching after structural/coverage
filtering is allowed. See
[the expanded workflow](terminal-sister-backgrounds-20260927.md#expanded-catalog-refresh-september-30).

### Expanded duplication backgrounds and reference choices — September 30

The complete expanded terminal-pair background support workflow retains all
140,719 native family trees and 3,154,373 pair records across the two guides.
Full source/readback checks cover 160,415 modeled candidates, 146,172 distinct
eligible model pairs, four annotation policies and all 3,400,908 architecture
support rows. The target universe remains all 283,409 modeled duplicate links;
no passing-screen-only or outcome-selected subset defines the matching pool.
All 5,059,122 candidate edges and both possible endpoint metadata mappings
passed independent checks. Fixed selection retains the original 54 scenarios,
zero-distance rule, deterministic ties and unmatched/reuse dispositions,
requiring 61,216,344 decisions. Full independent readback subsequently passed all decisions, including 4,250,692 selections and 56,965,652 unmatched outcomes.
Controls will not be rematched after structural or coverage filtering.

Separately, expanded provisional sister references were independently
reconstructed from each native tree for every target. All nearest modeled
ties and both duplicate sides yield 121,490 ledger rows; 86,440 native
gene/model joins and every exported catalog descriptor passed independent
readback. Cross-guide comparison retains the full 142,107 gene-pair union:
28,950 eligible targets share the same nearest-reference gene set and 66 have
disjoint sets. Extant references do not reconstruct ancestral states; their
availability-dependent selection requires sequence-locked sensitivity and
biological orthology assessment before asymmetry inference.

A full pair-union catalog source inventory compares old/new primary,
reference and background inventories and includes the already measured
expanded primary queue. An independent SQL join checks every pair and source
membership, including numeric model versions. Matching coordinate/sequence
checksums and lengths identify planning candidates only. Actual reuse requires
raw coordinate and materialized input/mask/order equality, executable/options
and checkpoint bindings, full numerical proofs and preserved exclusion flags.
Full additional-background raw-coordinate validation is active. Expanded
measurements, balance, prediction uncertainty, phylogenetic dependence and
inferential calibration remain incomplete. See the
[background workflow](terminal-sister-backgrounds-20260927.md#expanded-support-and-candidate-graph-completed-september-30)
and [reference workflow](duplication-sister-references-20260926.md#expanded-catalog-reference-pipeline-completed-september-30).

Full distinct-pair source inventory subsequently passed independent SQL checking:
336,292 current pairs, 236,650 matching existing catalog source signatures and
99,642 new to those collections. Matching-source candidates still require full
input/result and numerical checks before reuse. Full additional-reference
raw-coordinate validation and independent reconstruction are running for
24,804 models; these CPU stages do not resume protein prediction or qualify
asymmetry inference. All eight scientific aims remain incomplete.

Full descriptive pre-measurement matching balance passed across all
432 strata and eight features, preserving original target baselines, reuse and
nonestimable cases. Supplemental background/reference input preparation and
independent PDB/residue-map checking are queued behind complete coordinate
readbacks. They retain 607,604/49,608 two-mask dispositions, including short
and rejected inputs. Role overlap is explicit; future reuse requires exact
source, input bytes, mask and position equality. No result reuse or native
alignment is authorized by catalog signatures alone. Resource plans and
full handoff/corruption fixtures are versioned in
`metadata/expanded_additional_alignment_inputs_queued_20260930.json` and
`metadata/expanded_additional_alignment_input_fixture_validation_20260930.json`.

Expanded pre-measurement balance independently checked all 4,250,692 selected
records and 3,456 feature summaries; original-target selection shifts remain
explicit. In the unchanged illustrative S45 policy scenario, mean/minimum
pLDDT shift by approximately 0.843/0.911 original-target SD despite small
within-match standardized mean differences. This descriptive ascertainment
contrast requires sampling/missingness treatment; it is not prediction
accuracy, a calibrated effect or proof of causality. Complete strata and
source-checked PNG/PDF/SVG outputs are documented in
[the expanded balance record](terminal-sister-backgrounds-20260927.md#expanded-pre-measurement-balance-independently-verified).

### Sequence-first sister-reference availability sensitivity — September 30

For every modeled duplicate target, closest nonfocal sister genes were chosen
from the fixed native gene tree before consulting structure availability.
All tied genes within 1e-12, lexical representatives, missing models and original
parent exclusions were retained alongside the available-reference design.
Independent leaf-interval subtraction/upward fsum reconstruction checked all
283,409 records and context/choice summaries. Among 29,155 profile and 29,130
MAFFT provisionally available references, 7,560/7,566 use a farther gene because
the nearest sequence-only ties lack structures (25.93%/25.97%). Another 14/15
have an unmodeled lexical representative but a modeled equally nearest tie.
These are dependent conditional design diagnostics, not biological effects
or evidence of a preferred gene tree. Available-reference and sequence-first
ledgers remain separate; background controls are not rematched. Reference
orthology and complete structural measurement/reuse, coverage/confidence/
common-residue tests remain required.

Full additional-reference raw-CIF checks and written-input checks subsequently
completed: 24,804 models/9,982,283 residues and 49,608 mask dispositions, with
1,023 pLDDT70 inputs too short for alignment. All source, coordinate, original
position and written PDB checks passed, preserving short masks and source
provenance. Inputs are not pair results; full native result reuse verification
and numerical/biological qualification remain unfinished. See
[the reference evidence](duplication-sister-references-20260926.md#sequence-first-reference-choice-completed-september-30).


The completed sequence-first cross-guide sensitivity retained the full
142,107 duplicate gene-pair union, including parent-ineligible and one-guide
contexts. An independent dataframe outer merge checked every saved field and
context/choice matrix. Among 69,683 shared targets with eligible sequence-sister
contexts, 69,535 nearest gene sets agree and 148 are disjoint. Both lexical
choices have frozen models for 21,451 targets; 21,411 share a model/version.
Unmodeled/unmodeled is not model agreement, shared models do not establish
orthology, and overlapping guide records are not independent events.

The expanded reference alignment workload uses the full 55,701-pair ledger.
All 24,947 pairs without matching existing catalog sources undergo both input
orders and the original full/pLDDT70 masks (up to 99,788 dispositions); 30,754
matching-source pairs remain pending actual input/result/numerical reuse
qualification. The current model/input partition and full source workload
passed actual-data preflight before launch. Native measurements are active;
full least-squares/native-text diagnostics and rank/rotation-curvature and
quaternion checks are queued. Complete result union, original-length coverage,
common-residue triads, domain/PAE controls and calibrated asymmetry inference
remain unfinished. This resource partition is unrelated to structural outcomes
and does not authorize reuse or remove existing numerical flags.


The full reference reuse qualification checks all 55,701 model pairs/two masks/
two orders, including the 24,947 currently unmeasured pairs as pending states.
Completed expanded-primary measurements are preferred for 365 pairs and old
reference measurements for 30,389 pairs, without outcome-based selection.
Actual raw/PDB bytes, residue/sequence/mask/source signatures, executable
options and every directed checkpoint and original numerical/geometry flag
are checked before importing a compatible result. The producer records all
123,016 compatible old states with 118,804 usable and 4,212 unavailable or
numerically excluded directions; five RMSD flags remain explicit. The complete
independent reader and hash/journal closure passed all exported states, 348,748
source/artifact bindings and both exact captured process completion journals.
Compatible states retain all numerical exclusions for downstream union. Original results
and exclusions remain intact, and absent results are not treated as zero.
The full result union, order/coverage/common-residue/domain/PAE/orthology and
calibrated evolutionary analyses remain unfinished.

### Complete nearest-reference native orthology membership — September 30

Native orthology queries preserve all 283,409 source target contexts under both
tree methods, both reference designs, every nearest tie, missing model and
parent exclusion. The full design contains 214,461 reference tie records and
428,922 duplicate-to-reference links. Unordered physical gene queries are
deduplicated without deleting logical context/design/tie/side records. Each
reference is queried against both duplicates in both globally audited native
reciprocal streams. Protein ordinals are source-line positions, with mappings
bound through full protein/family identity audits and resolved-tree readbacks.
Compiled full-stream merge and independent fixed-record binary searches
provide separate membership checks. The reader independently reconstructs
every context, gene/model/version, lexical choice, eligibility and membership
summary, with bracketing for absence. Empty reference sets remain unqueried;
positive membership cannot remove a parent-context exclusion. Original fields
and both-guide coorthology states are retained. All 166,829 unique physical
queries and 428,922 logical links passed independent reconstruction; closure
checked 103 source/artifact hashes and both exact process journals. This is
native assignment provenance rather than independently validated orthology,
ancestral reconstruction or a structural duplication effect. Small-family
supplements are outside these native streams; all subsequent structural and
phylogenetic/calibration controls remain required.

Premeasurement cohort attrition was tabulated for every source context/design
and independently reconstructed using a uniquely keyed SQL table rather than
the producer's counter logic. All 566,818 rows, 48 nested eligibility/model/
native-assignment policies and 76 lexical context-state cells passed. Missing
genes remain unqueried, missing models are separate from native membership,
and any/all tied-reference sensitivities are reported without altering lexical
choices or parent exclusions. Parent-eligible lexical modeled references with
assignments to both duplicates under both native guides number 28,569/28,551
for availability and 21,169/21,150 for sequence-first designs. These overlapping
conditional records do not establish independent observations, quality-qualified
structural sample sizes or biological effects. Summary closure checked all 23
source/artifact bindings and two additional exact completion journals.

### Full reference measurement union — September 30

The full 55,701-pair/two-mask/two-order design retains 222,804 dispositions.
The join of 123,016 qualified old states and 99,788 new native states requires
complete new-native numerical/quaternion verification and all four original
process journals. Original checkpoint paths/checksums, native order,
ordered model/version endpoints, native metrics and complete numerical/
geometry records remain explicit. Different old native directions map by
actual ordered endpoints; no metric, flag or original checkpoint is changed.
All missing/error/short/degenerate/RMSD dispositions remain in the union, with
no zero imputation, optimization retry or outcome-based source substitution.
An independent reader reconstructs every exported field and source/status/
exclusion count before full byte/journal closure. Complete byte indexes stay
outside Git with versioned locations/checksums. These stages are queued behind
the live native measurements, not completed production results. Numerical
union does not establish original-length coverage, common-residue/sequence-
locked correspondence, domain/PAE/prediction accuracy, biological orthology,
accepted phylogenetic uncertainty or calibrated asymmetry.

### Complete reference input-order and original-length coverage — September 30

The complete reference design has 55,701 model pairs, two confidence masks
and two native input orders (222,804 dispositions). Each directed numerical
record is normalized to the current ledger's endpoint models and versions
using actual ordered endpoints; its original checkpoint order remains
unchanged and traceable. Native errors, unavailable masks and numerical
exclusions retain their original flags and raw aligned lengths, with blank
usable metrics. Reversed old source orders are not interpreted as reversed
current endpoint identities. Retained-input coverage and original-protein
coverage are separately named. Original full-protein lengths are checked
against frozen model descriptors and both complete written-input manifests,
including full sequence hashes and retained residue positions.

For each pair/mask, both numerically usable orders must meet each fixed screen:
30 or 50 aligned residues and at least 50, 70 or 90 percent of both original
protein lengths. Screens use exact integer fraction comparisons; the
independent checker reconstructs all source/metric fields and uses decimal
ceiling thresholds. All 111,402 pair/mask rows and 668,412 screen decisions
remain in the denominator, including excluded cases. Both-mask intersections
and absolute differences between usable orders are sensitivity summaries,
not independent replicates or uncertainty intervals.

Implementation and software fixtures are complete. Actual static input
preflight passed all 83,207 reference models and 602,972 frozen input states;
its 60 source/artifact bindings and original process journal are closed. Full
production, independent readback and journal closure are queued behind the
closed measurement union. Event/context linkage, residue correspondence and
common-triad/sequence-locked/domain/PAE/prediction, phylogenetic and calibrated
inference controls remain outstanding. Coverage qualification alone does not
accept a reference's biological role or establish structural asymmetry.

### Complete native context-to-measurement design — September 30

Every one of the 283,409 closed native target contexts was joined to the
frozen duplicate gene/model/version queue. Gene identity, rather than table
a/b position, determines each model role. All native context/source fields,
parent eligibility, original reference designs, tied genes, lexical choices
and both-guide native assignments are retained. For each of the 214,461
reference ties, both duplicate sides are mapped to exact current reference
pair endpoints/versions, or retained as missing-reference-model, identical-
model or outside-design work states. All 121,490 original availability-side
ledger links are accounted for. Matching physical pairs may be shared among
events, designs or reference gene aliases, without deleting logical links.

The full 428,922 logical-side export and 566,818 context/design denominator
were independently reconstructed with uniquely keyed SQL queue/pair/ledger
joins. Every source/native field, gene/model/version role, pair hash, focal
endpoint, missing/identical/outside-design state, lexical flag and parent
exclusion passed; all source ordinals and both full context/availability
universes were exhausted. Completion verified 40 source/artifact hashes and
both exact original process journals. Work-state counts use explicit counting
units: no-reference cells count contexts, while other cells count logical
duplicate/reference sides. These units are not combined.

This completes the fixed source-only bridge, not structural cohort acceptance.
Original parent eligibility must still gate every subsequent coverage/native-
assignment policy, including when an excluded context maps to a physical pair
measured for another event. No modeled tie replaces an unmodeled lexical
choice, no missing comparison becomes zero, and an identical model is not an
independent structural replicate. Coverage, common-residue and sequence-locked
comparisons, domain/PAE/prediction uncertainty, biological orthology/phylogeny
and calibrated inference remain required.

### Full context coverage and native-assignment policies queued — September 30

The complete fixed context design is joined to the source/journal-closed
reference pair/mask coverage matrix. All 283,409 context records and every
source/native field remain in the full export. Both reference designs, both
guides, all 214,461 tied-reference records/428,922 duplicate-reference sides,
parent exclusions, missing structures and identical/outside-design states
are retained. Measured pair flags use both numerically usable native orders
and the original full-protein denominator; unmeasured order status/numerical
fields are null and their screen exclusions are categorical. Original
measured numerical flags are copied unchanged from the verified matrix.

For each design/mask/screen, both duplicate-to-reference sides must pass
coverage. Original parent eligibility gates every context policy. Five
policies report lexical reference coverage alone, lexical coverage with
assignment to both duplicates under the source guide, lexical coverage with
assignments under both guides, any tied reference meeting coverage and both-
guide assignment, and all tied references meeting those requirements. Empty
reference sets fail every policy. Any/all-tie sensitivities are separately
reported; they never replace the source lexical gene. Physically measured
pairs do not override excluded parents. These diagnostics do not yet qualify
primary duplicate AB or common-triad independence/residue correspondence.

Both masks and six original screens represent 6,801,816 context/screen states
and 34,009,080 contextual policy decisions, with 5,147,064 logical-side/screen
decisions. These are stored in the source-preserving nested compressed
export, not interpreted as independent observations. The 240-cell summary
retains source-context, parent-eligible, parent-eligible lexical-gene and
parent-eligible lexical-model denominators. An independent checker rebuilds
each side state using SQL pair lookups, derives all contextual flags with
SQL boolean/count aggregates, and checks all source ordinals, dimensions,
summary cells and denominators.

Code and full software fixtures passed, including excluded parents with
measured pairs, a missing lexical model with modeled alternate, empty all-
tie sets, native guide disagreement and numerical exclusions. Eight rehashed
false exports were rejected. Actual context source/40 hashes/current inventory
and primary-queue bindings passed preflight. Full production, SQL readback
and both original process journals are queued behind closed coverage, not
completed results. Primary AB/common-triad/sequence-locked/domain/PAE/
prediction, biological orthology/phylogenetic and calibration checks remain
required before effects can be tested or interpreted.

### Full duplicate/reference three-edge work design — September 30

The complete closed context/measurement design was extended without dropping
any of its 283,409 source contexts, 214,461 reference ties or 428,922 logical
duplicate/reference sides. Original parent eligibility, gene orientation,
lexical choices, model availability and both-guide native assignments remain
unchanged. Primary AB and full-reference AR/BR pair catalogs determine actual
current endpoint orders and desired A→B, reference→A and reference→B directions.
Unmeasured/identical endpoint-order fields are null. Modeled references define
ordered versioned A/B/reference triples; repeated physical triples are
deduplicated without merging their logical source occurrences.

Source-only correspondence work requires an eligible original parent, queued
primary AB, all three pair-catalog memberships, three distinct versioned models
and three distinct model IDs. Version differences alone do not establish
independent proteins. Own-guide and both-guide native assignment sensitivities
are separately reported. This source gate does not assess native measurement
coverage or biological orthology. Independent SQL catalog queries reconstruct
every exported field/edge/direction and full logical/physical link counts;
source ordinals, identities, both designs and all empty/missing/excluded states
are checked. Full independent readback passed and exact two-journal closure
verified 59 source/artifact hashes.

There are 31,235 unique ordered model triples; 27,056 have at least one
source-ready context/reference occurrence. These define 432,896 future
correspondence states under two masks and all eight AB/AR/BR native order
combinations. No residue-map or common-core fit is claimed here. The next
stage must use actual native source directions and original residue positions,
retain all numerical/error exclusions, compare reference-common versus
AB-cycle-consistent cores, fit identical residue triples and evaluate original
full-protein coverage. Sequence-locked/domain/PAE/orientation, prediction-error,
phylogenetic and calibrated inference controls remain required. Exact count
units and provenance are in [the full workflow](full-reference-triad-work-design-20260930.md).


### Full original-residue duplicate/reference correspondence

The completed source-ready grid includes 27,056 ordered physical triples,
both masks and all eight AB/AR/BR native orders, without coverage preselection.
Actual raw checkpoint inputs, directions, source orders and checksums determine
pair mappings on original protein positions. AR/BR intersections and a separate
AB-cycle-consistent subset are retained, including failed, numerically excluded
and short states. Both-order pair-screen flags retain original full-protein
coverage denominators. An independent iterator/relation reader reconstructs all
432,896 states and every classification. Exact two-journal closure binds 321,723
hashes. No shared-coordinate fit or evolutionary effect is established by these
checks. Input sources, complete contextual coverage and limitations are recorded
in [the full correspondence workflow](full-triad-common-residues-20260930.md).


### Full identical-residue geometric comparisons

The scheduled 865,792 fit dispositions retain both masks, eight edge orders and
reference-common versus cycle-consistent cores. Proper SVD rigid fits for AB,
AR and BR use identical original-position triples; an independent quaternion
reader reconstructs every distance, orientation-preserving contrast, confidence,
identity and uniqueness classification. Short/source-excluded fits have blank
metrics and degenerate computed fits stay excluded. Six rational original-length
core screens and inherited both-order three-edge screens are reported separately
and jointly. Full source/journal closure gates production. The exhaustive work
preflight measured 77,805 distinct eligible cores and 233,415 pair fits per
implementation. The complete numerical geometry subsequently passed independent
reconstruction; [production and resource details](full-triad-same-residue-geometry-20261001.md).
Calibrated asymmetry is not yet claimed.


### Full crossed PMSF topology, consensus and character diagnostics

All four fixed alignment/guide combinations retain both ML and consensus trees.
Six cross-run pairs per tree type and four within-run pairs exhaust every split,
missing cell and incompatible split pair on the same 526 tips. ML support
criteria use SH-aLRT80/UFB95; consensus uses UFB95 with SH-aLRT unavailable.
DendroPy independently reconstructs raw nodes/lengths, compatibility and quartet
witnesses. Comparisons retain all support alternatives rather than selecting a
preferred run. Full taxa-character coverage and nearest declared outgroup-set
boundaries retain every taxon and tied minimum. Counts use independent complete
FASTA parsers and topology libraries. These diagnostic distances are not rates,
biological discordance tests or proof of a root. Complete results, warning
interpretation and numerical/journal evidence are in
[the species-tree sensitivity workflow](pmsf-four-run-sensitivity-20261001.md).

The full same-residue geometry subsequently passed all 865,792 dispositions and
independent quaternion checks with maximum difference 5.107e-14 Å. Core and
inherited three-edge six-screen passes remain separate; logical events and
order/mask robustness are not inferred from aggregated repeated state counts.


### Full structural order robustness and original context linkage

All 865,792 geometric dispositions are aggregated into 108,224 physical
triad/mask/core groups, retaining every eight-order status, residue hash,
exclusion and combined/core/inherited screen bitmap. Eighteen numeric fields
retain available-order counts and ranges; missing metrics remain null. Strict
all-order sign categories keep original model roles and are descriptive, not
calibrated evolutionary directions. Independent complete SQLite reconstruction
checks every source row/group/summary; maximum numeric difference is 2.842e-14.

The full context projection embeds unchanged original work records and links
all 283,409 contexts/214,461 tied references to normalized physical results.
Original source/parent/model-identity and native-guide gates remain mandatory,
including measured overlaps in otherwise excluded contexts. Both confidence
masks and their exact order-bit intersection remain separate for each core.
Five fixed lexical/native/any/all-tie policies retain missing references and
empty ties without reselection. Full model-role/source-gate and SQL policy
reconstruction passed every contextual output and summary, with 464,637 source
bindings and both original journals. Counts use
physical triples and repeated logical quality screens explicitly; they are
not independent events. [Full workflow, evidence and limitations](full-triad-order-context-robustness-20261001.md).


### Full expanded background measurement workflow

The 146,172-pair background design retains 305,434 active models and both
confidence masks. Exact physical input collection order is declared in advance
across primary, modern reference, background and legacy reference collections;
every global overlap requires identical full semantic fields and actual written
bytes. All origins and original residue positions remain traceable. A full
independent source/SQLite readback and two original completion journals must
close this union before native work. All 29,080 legacy written-input states
have passed readback and two-journal closure; production union is still running.

The queued native stage measures every 74,712 new pair under both masks/orders
(298,848 dispositions). The 71,460 catalog-matched pairs retain pending reuse
status until exact old-input/result/numeric/geometry checks succeed. A complete
checkpoint assessment preserves nonalignment and numerical exclusions and
reconstructs least-squares distances, sequence identities, confidence, coverage,
coordinate rank and proper-rotation curvature. Separate vectorized mapping,
quaternion distance and LAPACK spectrum checks audit all states. Three exact
original native/assessment/reader journals and full source/artifact hashes are
required for new measurement closure. Software checks passed on actual native
synthetic coordinates and rejected 12 rehashed false outputs; this is not
production completion or a biological pilot. TM-scores bind native output
without independent optimization. No biological effect, rate or calibrated
selection inference follows from numerical eligibility alone. [Full workflow
and prelaunch resources](terminal-sister-backgrounds-20260927.md#full-four-collection-input-union-and-new-measurements-october-1).

Subsequently, full input-union source/SQL readback and both original journals
passed, binding 964,012 source/artifact hashes. All 610,868 active input states,
global overlaps and 146,172 work pairs are closed. The native runner has started
input checks; production native numerical/geometry completion remains pending.
Actual old-result qualification is launched over all 584,688 full-design states.
Reference-before-background preference is fixed, with 11/71,449 source pairs;
new, incompatible and nonalignment states remain explicit. Source coordinate,
original residue/sequence/length/mask/status/PDB bytes and native settings must
match. Every compatible original numeric row is reconstructed on the current
coordinates; separate quaternion and alternate-SVD checks qualify the complete
reuse exports and old geometry, with original numerical exclusions retained.
Eight original source journals passed separate verification; final reuse
closure requires both new plus eight original journals. Actual-native synthetic
software checks rejected 12 rehashed exports, preserving reversed orders,
incompatibilities and RMSD/nonunique/unavailable states. Production reuse,
new/reused result union and calibrated evolutionary inference remain pending.


The complete background result union is now queued over all 146,172 physical
pairs/584,688 mask-order states. It requires independently closed new-native
(three original journals) and actual old-result reuse (ten original journals)
sources. The reuse producer has finished, reporting 285,840 retained old states
and 298,848 new states pending without incompatible inputs; full independent
validation remains active. Source ownership and original alignment direction
are fixed before numerical outcome. The union retains original numeric and
geometry rows beside typed fields, every native nullable/error/timeout record,
model roles and all exclusions. A separate full SQL reader reconstructs each
field and count from the closed sources; final closure rechecks all bytes and
both original union journals. Native optimization/geometry are not rerun at
this merge. Software fixtures rejected 15 rehashed false exports and blocked
incompatible sources; synthetic prior measurement/journal contracts do not
constitute production evidence. [Scope, resources and reproduction](terminal-sister-backgrounds-20260927.md#full-background-measurement-union-queued-october-1).
Production completion, coverage/prediction/phylogenetic controls and calibrated
evolutionary effects remain pending.


The full background original-protein coverage stage is queued behind successful
whole-measurement union closure. All 146,172 physical pairs/292,344 pair-mask
rows/584,688 original states enter the same six target thresholds, yielding
1,754,064 decisions. Both alignment directions must be numerically usable;
30/50-residue and 50/70/90-percent cutoffs use both original full protein lengths
under either confidence mask. All original statuses, flags, sources, directions,
nullable values and overlapping exclusions remain explicit. Full 305,434-model
input inventory and 292,326 background endpoint identities are counted
separately. Independent SQL-original-state/decimal-ceiling reconstruction checks
every row, intersection and 36 target/background physical count rows; completion
requires all bytes and two original coverage journals. Software fixtures rejected
13 rehashed false exports and cover unused inventory models and exact cutoff
boundaries; synthetic prior proofs are not production evidence. All fixed
selected/unmatched scenario integration and post-screen balance remain downstream,
without rematching. [Scope and resources](terminal-sister-backgrounds-20260927.md#complete-background-coverage-screen-queued-october-1).
Independent old-result readback has passed all 584,688 design states; its full
ten-journal closure and new-native measurement/union completion remain pending.


## Full frozen matching coverage attrition — October 1

Before integrating responses, we checked all 283,409 duplication-target and
318,037 background graph nodes against their original versioned model
catalogues and the appropriate complete physical comparison tables. Target
models were not required to occur in the background-only active inventory.
All roles, model IDs/versions, sequence and coordinate hashes, original lengths
and catalogue confidence values were checked, including same-model nodes.
The graph/covariate four-journal closure was bound separately from the
selected-control two-journal closure, with explicit producer/readback and
fixed-selection lineage checks.

The queued integration keeps the original 54-scenario/four-policy/two-guide
matching design fixed. Original selection identities, endpoint order, scores
and ties, plus all matched/unmatched lists, are retained. Every physical
comparison receives the unchanged six target coverage screens (30/50 aligned
residues and 50/70/90% of both original full proteins), separately for full
and pLDDT70 masks and their intersection. Same-model comparisons receive an
explicit no-alignment exclusion, not a zero response. All selected failures
and unmatched decisions remain in their original denominators.

An independent reader compares all 4,250,692 selected source/export rows and
1,133,636 target-policy rows, independently projects flags by endpoint
identity and reconstructs every one of 7,776 attrition cells in SQLite,
including empty strata and all joint/target-only/control-only/neither-pass
categories. Full closed input hashes, all output hashes and both original
completion journals are required before production acceptance. Software
fixtures cover all 54 scenarios and reject corrupted exports/provenance;
their prior physics and journal contracts are synthetic. Production results
remain pending. Counts are dependent design records; post-screen covariate
balance and calibrated phylogenetic duplication effects are separate analyses.


## Full post-screen balance and control reuse — October 1

The queued post-screen analysis retains every original fixed-match scenario
and all closed source denominators. Eight permutation-invariant sequence,
length and confidence features are summarized among jointly eligible target/
control records. We distinguish shifts relative to all original modeled
targets, original metadata-matched targets and all target-screen-eligible
records (including unmatched targets), preserving separate finite-feature
counts. Zero-distance log exclusions, empty cells, exact constant variance
and insufficient pairs/baselines remain explicit. No epsilon, outcome-based
rematching or balance-pass threshold is introduced.

Full representation summaries count target taxa, families, genes and physical
pairs, and control taxa, genes, nodes and physical pairs. A retained reuse
ledger reports node and physical-pair multiplicity separately, with reciprocal
weights assigning one total unit per control under each definition. Weight
sums and Kish concentration are descriptive diagnostics, not phylogenetic
effective sample sizes or accepted inferential weights. Independent scalar/
pandas reconstruction checks all 62,208 feature rows, 7,776 coverage cells
and all reuse records, including empty and nonestimable cases, followed by
full source/artifact/two-original-journal closure. Software fixtures use
synthetic prior matching/physics/journals; production results are pending.
Biological duplication effects require subsequent dependence, predictor,
structural-domain/PAE/orientation, missingness and phylogenetic calibration.


## Native taxon species-tree sensitivities — October 1

Four taxon policies are applied separately to the original profile and MAFFT
concatenations: exclusion of sparse boundary taxon F1243177; canonical occupancy
of at least 10% in each alignment; exclusion of two assembly-linked curated
hybrids; and exclusion of hybrids plus 21 incomplete labels. Exact retained
sequences, all 125 markers, columns and partition/site mappings remain unchanged.
Independent full character reconstruction verifies each subset and all membership/
lineage counts. These cohorts do not replace the full baseline or establish
accepted unique-species counts. Their loss of outgroups and fine taxonomic
coverage is retained for interpretation.

All four policies cross both matrices with both conditioning guide alignments
in 16 supported LG+C20+F+G4 PMSF refits. Four byte-identical existing homogeneous
guide outputs are checked as frozen conditional inputs; missing guides are
freshly inferred from the complete corresponding sensitivity data. Every PMSF
profile and supported tree is refitted with 1,000 SH-aLRT/1,000 UFBoot replicates
and BNNI, then receives full profile/report/raw-tree/empirical-bootstrap/NEXUS
readback on its exact taxon universe. Runs are serial with a memory gate;
complete collection reconstruction, common-tip comparisons, root/marker/model/
identity sensitivity and original-journal closure remain required. This is
ongoing sensitivity inference, not a pruned-tree substitute or final phylogeny.


The native taxon-sensitivity outputs will additionally undergo independent
DendroPy tree parsing and Decimal profile reconstruction. For each of the
16 crossed runs, the reader enumerates all raw bootstrap tip/edge grids,
reconstructs empirical split frequencies and checks the full NEXUS universe,
raw/report tree agreement and all support fields. It uses the actual retained
manifest roles for each ML/consensus boundary diagnostic. Original batch
completion and independent reader completion must be demonstrated by their
captured PID/create/command journals and full source hashes. Software
regression against a complete existing baseline and synthetic I/O/corruption
cases is complete; the actual sensitivity readback is still pending. These
checks do not independently recompute SH-aLRT tests or likelihood optima and
do not establish model adequacy, a root or an accepted species framework.


For comparisons between native taxon-sensitivity refits and full-cohort
inference, baseline trees are projected onto each exact retained cohort.
Empirical ultrafast-bootstrap frequencies are recalculated by counting each
projected split once in each original replicate, even when multiple original
edges merge. Original edge lengths are summed along retained paths; these
values remain conditional on full-cohort inference. SH-aLRT support is not
inherited. Independent taxon deletion and path suppression in DendroPy,
including a unit-edge reconstruction of path component counts, checks the
complete projected edge/frequency/role tables. Production and original
journal closure are still pending. This diagnostic supplies a matched-taxon
reference and does not replace actual subset inference or establish a root,
model adequacy, dates or calibrated sequence–structure effects.


Projected baseline topological sensitivity will be compared within each exact
retained cohort across all crossed alignment/guide runs and ML/consensus
views. Every incompatible pair of unique splits is retained with a quartet
witness, and every union split receives a presence/absence row for each
view. Recomputed projected UFB >=95 is a descriptive support screen; SH-aLRT
is unavailable. Independent raw-tree pruning, DendroPy RF and bipartition
compatibility checks reconstruct all comparisons and conflicts. Production
and journal closure remain pending, and native-versus-reference comparisons
will remain necessary. Reference-only comparisons conditional on the original
full-cohort fit do not establish native subset robustness or model adequacy.


### Full gene concordance across candidate species views (October 2 UTC)

All 125 profile marker trees and all 125 MAFFT marker trees were compared with
each of the eight closed original PMSF ML/consensus views using IQ-TREE 3.0.1
`--gcf --cf-verbose`: 16 runs, 8,368 branch summaries and 1,046,000 branch–gene
cells. Decisiveness requires coverage of all four incident reference clades.
Concordant, two NNI alternatives, residual discordance and uninformative
missing-clade states remain separate, with the actual decisive denominator.
Reference internal labels are removed before native reorientation; original
SH-aLRT/bootstrap support is joined by canonical split.

Independent DendroPy reconstruction checks every restricted gene split,
incident-clade count, native verbose/aggregate/NEXUS value, branch ID, factor
label and branch length. Native NNI orientation must have one consistent
mapping across genes for each branch. Zero decisive genes retain unavailable
factors and all-NA native statistics. Numerical checks respect two-decimal
statistical percentages, three-significant-digit Newick factors and six-
significant-digit statistical branch lengths. The initial reader inventory
failure is retained; new v2 output checks every NEXUS annotation. Full v2
readback and two-original-journal closure passed all 1,046,000 cells with
1,217 bindings. Descriptive support joins retain all 8,368 branch rows,
including ten unavailable factors; all counts/factors and 16 medians were
checked against the closed reader ledgers and the figure was inspected.

See [complete workflow and interpretation](species-gene-concordance-20261002.md)
and [Minh, Hahn and Lanfear 2020](https://doi.org/10.1093/molbev/msaa106). Raw
gCF is conditional on inferred marker/reference topologies and taxon coverage;
it does not distinguish biological discordance from estimation error or
qualify a root, reconciliation, dating or structural evolutionary effect.


The retained-reference projection and comparison stages subsequently passed
full independent reconstruction and two-original-journal closures: all
16,000 projected bootstrap states, 32 views, 33,088 edges and 19,289 empirical
split-frequency rows (207 bindings), followed by all 64 comparisons, 17,744
presence cells and 1,898 overlapping incompatible split pairs/quartet witnesses
(229 bindings). Closed small receipts locate the complete source archives:
[projection](../metadata/retained_taxon_projections_completed_20261001.json)
and [comparison](../metadata/retained_tree_comparisons_completed_20261001.json).
Original path sums and missing projected SH-aLRT remain explicit; full actual
native subset refits and native-versus-reference comparisons are still needed.


### Full coalescent species-tree sensitivities (October 2 UTC)

ASTRAL-III 5.7.8 completed every combination of profile/MAFFT genes, five
full/closed taxon cohorts and three support settings: uncontracted, SH-aLRT
below 10 contracted, and SH-aLRT below 80 contracted. All 125 markers remain
in every input. Missing support is recorded separately and contracted in
filtered inputs; exactly 10/80 is retained. Contraction precedes pruning.
These are SH-aLRT sensitivities, distinct from published bootstrap thresholds.
The primary set keeps all 526 entries and 25 outgroups; species identity
qualification remains open. Full source census passed all 250 trees, ten
cohort/alignment unions and 1,305 bindings; every retained marker has at
least 402 taxa.

All 3,750 input states, 1,781,685 support decisions and 1,545,528 retained
splits passed independent DendroPy reconstruction. Closure verified 1,408
bindings and both original journals. Installed jar and three dependencies
match the official 5.7.8 archive. Native resolved/missing/polytomous contracts
passed. All native commands, child process identities, return codes and
annotated tree/log hashes are retained outside Git.

Independent Python/Numba readback passed every one of 15,510 branches and
1,938,750 branch–gene states. A colored postorder dynamic program computes
exact local quartet counts; gene-specific normalization preserves missing
clades and unresolved evidence. Separate node-component intersections
compute the exact global matching score, and unresolved quartet centers
give the independent global denominator. Beta-tail integration and the
native default lambda = 0.5 recompute local posterior support and MAP
coalescent lengths. NNI alternatives are compared as an unordered pair.
All 18,720 color/tree cases, five native numerical contracts and 1,264 global
tree-pair cases passed; 120 altered native values were rejected.

The first actual full candidate passed all 523 branches and 65,375 local
states, its exact global numerator/denominator and one-journal closure with
1,434 bindings. This certifies one candidate's numerical calculations; the
full 30-case readback subsequently passed all local/global states and
two-journal closure with 1,757 source/artifact bindings. All 315 comparisons
against concatenated references also passed. Full readback resources were
estimated before launch: two CPU/16 GiB/no swap/one BLAS thread, no GPU or charges.

[Workflow, complete inputs and current evidence](species-coalescent-sensitivities-20261002.md).
Candidate estimation remains conditional on inferred genes, taxa/support
settings and MSC/locality assumptions. Numerical agreement does not prove
search optimality or biological model adequacy. Native rooting is arbitrary;
coalescent lengths do not estimate substitutions, time or structural change.

### Matched-taxon coalescent/reference tree comparisons (October 2 UTC)

Thirty numerically closed coalescent candidates are compared with the 40
original/projected concatenated ML/consensus views under five exact cohorts.
Every six-candidate/eight-reference cross contributes 48 pairs; all 15
coalescent candidate pairs per cohort are retained, total 315 comparisons.
Every union split-presence cell, RF count, incompatible unique-split pair,
canonical quartet witness and all 70 ingroup/outgroup boundaries are exported.
An independent DendroPy raw-tree pruning/RF/bipartition reader reconstructs
every export before full hashes and two original journals permit closure.

Local posterior 0.95, original ML SH80/UFB95, consensus UFB95 and projected
UFB95 with unavailable SH-aLRT remain distinct descriptive criteria. Effective
genes and branch units are retained. Coalescent MAP lengths are not compared
numerically with amino-acid substitution or original projected path lengths.
Role-boundary presence is unrooted and does not assign biological polarity.

The full synthetic 70-view/315-pair grid passed with 13 false exports rejected.
The actual first full candidate versus eight original references passed
5,814 presence cells and 1,988 overlapping incompatible pairs, with both
original journals and 1,475 bindings. Pairwise shared splits span 431–438/523;
417 splits are shared across all nine views. Comparisons condition on
inferred genes and candidate trees; overlaps are not independent events
or evidence selecting a biological cause or preferred model. Full 315-pair
production/readback/closure subsequently passed all 48,412 presence cells,
49,847 overlapping conflicts and 70 role-boundary rows (64 present), with
1,805 bindings and both original journals. The primary cohort shares 397
splits across all 14 views; the four retained cohorts share 410/411/399/384
respectively. Across all cohorts, cross-reference RF ranges from 134 to 184,
and within-coalescent RF from 14 to 138. These are conditional descriptive
measurements rather than qualified structural evolutionary branches.
[Workflow and limitations](coalescent-reference-comparisons-20261002.md).

### Complete candidate-tree sensitivity figure (October 2 UTC)

After full source comparison closure, all 240 coalescent/reference and
75 coalescent alignment/support pairs are displayed in ten cohort panels
on two PDF pages with a common color scale. Normalized RF uses twice the
retained taxon count minus six; every compared pair uses identical taxa.
Each pair remains linked to its source cohort/view identities and precise
distance. An independent reader checks every cell placement and all printed
three-decimal PDF distances, count identities, roles and unavailable triangular
cells. Full hashes and both original producer/reader completion journals
precede a separate visual inspection of both actual pages. The full synthetic
315-cell/two-page contract passed, rejecting 12 altered exports. Actual
production rendering and independent readback passed all 315 placements and
printed PDF values, followed by 1,833-binding/two-journal closure. Both actual
pages were visually inspected and published as a standalone PDF and PNGs;
no fixture is presented as biological evidence. Serial two CPU/8 GiB/no swap/one BLAS
resources were estimated before launch, without GPU or charges. These
descriptive point-tree differences inform later uncertainty analysis and
do not qualify a biological root, preferred model or structural effect.

### Native effective-N semantics and dependency recovery (October 2 UTC)

All 30 native tree outputs are complete and their 1,571-binding/original-journal
inventory closure is verified. The initial full numerical reader failed on
125.0 native effective N versus 124.99995265575164 fractional resolved evidence.
ASTRAL's implementation retains the count of genes having all four incident
clades when the resolved total differs by at most 0.001; otherwise it uses that
fractional total. The v3 reader keeps both totals, records this disposition,
and uses the independently verified native effective denominator for quartet
fractions, local posterior and MAP length. The exact quartet DP and 2e-8
absolute/relative comparison tolerances are unchanged. Three actual installed
5.7.8 boundary contracts and all five earlier native contracts passed;
120 altered values of 0.0001 were rejected. The previously failing complete
candidate passed all 65,375 states and exact global score; the full 30-case
batch and two-journal closure subsequently passed. Original failures/partial
exports are preserved, and new v2 plans/output directories carry the corrected full numerical/comparison/
figure workflow. Successful dependency handoffs require original invocation
completion/resource journals; collected unit defaults do not prove success.
[Native semantics and evidence](species-coalescent-sensitivities-20261002.md#native-completion-and-effective-n-correction-october-2-utc).

### Actual native taxon-subset tree comparisons (queued October 2 UTC)

The 16 actual subset PMSF fits will contribute 32 native ML/consensus views
under four previously fixed retained cohorts. Every native view is compared
with all matching projected concatenated references and coalescent candidates,
and with every other native view: 88 total views and 560 pairs, consisting of
192 native/coalescent, 256 native/projected and 112 native/native comparisons.
Every split-presence cell, incompatible pair, quartet witness and all 88 role
boundaries are retained. Native consensus multifurcations remain unresolved.
Observed split counts and missing internal slots are reported alongside
RF/[2*(taxa-3)] and RF divided by the observed split total; the latter is
unavailable for a zero/zero denominator. Independent DendroPy raw-tree
parsing/pruning/RF/compatibility/metric readback precedes full hash and both
original producer/reader journal closure. Native SH80/UFB95, consensus and
projected UFB95 with unavailable SH, and coalescent PP95 remain distinct
screens and branch units remain explicit. Full 88-view/560-pair software
contracts passed, including four zero denominators and 15 altered exports
rejected. Production waits for complete original native inference/collection;
these software contracts do not constitute biological results. Resources were
estimated before launch at serial two CPU/16 GiB/no swap/one BLAS thread,
without GPU or charges. [Workflow and limitations](native-subset-tree-comparisons-20261002.md).

### Full sequence-derived duplicate/reference correspondence (October 2 UTC)

The original full triad design was independently inventoried using full
sequence strings and exact model/version identities before confidence filtering.
All 31,235 original physical triples remain represented; 27,056 source-ready
triples are scheduled and 4,179 retain explicit unscheduled dispositions.
The catalog is source-bound to all 283,409 contexts, 214,461 ties and 428,922
logical sides. Every A/B/reference role assignment is preserved under sorted
model indices. Full independent source/sequence/role readback and two-journal
closure passed 321,761 bindings. Readiness does not establish biological
duplication or orthology.

Both installed MAFFT 7.525 (`--amino --anysymbol --thread 1 --auto`) and FAMSA
2.5.2-2598410 (`-t 1`, default single-input settings) are running for all six
input permutations of every set: 324,672 native dispositions. Common fully
observed columns yield exact original residue triples; no structural distance
selects an alignment or input order. Exact input/output/error bytes, checksums,
commands, PID/create observations and every native timeout/error/malformed
disposition are retained in transactional SQLite checkpoints. The separate
reader uses Bio.SeqIO and independent residue counters to reconstruct all raw
MSAs and their coordinate correspondences. Full native source/hash/readback
and two-original-journal closure remain pending. Installed-tool software
contracts passed all 48 indel/repeat/identical/nonstandard-letter states under
the complete two-method/six-order grid, rejecting 15 altered exports; exact
interrupted-task input replay also passed.

Resources preceded launch: two one-thread CPU workers, 32 GiB RAM, no swap,
32 GiB output allowance, 100 GiB disk reserve and no GPU or charges. Full
preflight, independent raw-MSA preflight readback, and their original-journal
closure are queued after complete native closure. The preflight measures
all 649,344 intended full/pLDDT70 fit dispositions and exact unique eligible
cores/PDB demand before coordinate fitting. Full SVD fitting, independent
raw-MSA/quaternion reconstruction and two-journal closure are also queued.
Both masks use identical three-way original residue matches, with joint
membership in authoritative confidence-mask position lists. PDB confidence
rounding never selects a residue; rounded-PDB confidence summaries and
authoritative-mask fractions are exported separately. Original-length and
retained-mask coverage, all six rational core screens, inherited both-order
three-pair exclusions and their intersection remain distinct. Short/degenerate/
source-rejected/native-failed states retain explicit dispositions and blank
uncomputed metrics. An exhaustive 432,896-view source check established
constancy of inherited both-order gates across structural-order views.

The complete synthetic 72-state/144-fit software grid passed real raw-MSA/PDB/
SQLite and SVD/quaternion algorithms, including nonconsecutive masks, a rounded
69.996-to-70.00 confidence boundary, unequal original lengths, duplicate-role
reversal, reflections, collinearity, short cores, native failures and inherited
exclusions. All 16 rehashed geometry exports and three changed preflight
summaries were rejected; interrupted recovery rechecked 24 committed rows.
Production sources are not replaced by these fixtures. Transactional fit
checkpoints are reconstructed and checked on restart. Full production preflight,
fits and independent readback remain pending. Prelaunch planning ranges are
uncalibrated and exclude dependency waits; no GPU or additional charges.
Comparison against structural correspondence, domain/orientation/PAE/predictor
controls and phylogenetic/calibrated inference remain required. These dependent
alternatives do not increase biological replication; alignment agreement does
not establish homology truth, prediction independence or ancestral polarity.
[Full workflow and current evidence](full-triad-sequence-correspondence-20261002.md).

### Full sequence versus structural correspondence sensitivity (October 2 UTC)

After complete source geometry closure, all 108,224 sequence method/mask/triad
groups are compared with both original structural core definitions. Each of
216,448 comparison groups retains all six-by-eight input-order alternatives:
10,389,504 paired states. The complete implementation, independent reader and
original-journal closure are queued; production results remain pending.
Every paired state retains original residue-triple intersection/union/Jaccard,
correspondence equality, both fit statuses and numerical uniqueness, strict
numerical sign agreement, all 18 shared numeric differences, and each source's
six combined screens with their intersection. Empty-union overlap and uncomputed
numeric differences remain unavailable. Different cores may contain different
residues or coverage; the differences quantify descriptive sensitivity rather
than independent biological effects.

Full sequence-order summaries retain source/native/payload/correspondence identity,
all 19 numeric ranges, core/inherited/combined screens and all six exclusions.
Joint summaries retain every 48-pair layout, range, missing state, sign count
and all/any-order bitmap. No favorable order or method is selected and alternatives
are not replicates. The independent reader rebuilds raw sequence correspondences
with Bio.SeqIO, merges sorted tuples for overlap, and reconstructs every pair,
range, mean, status and bitmap using SQL. Absolute numerical tolerance is 1e−9.
Full source/hash and both original producer/reader journal closure are required.
Deterministic per-triple checkpoints are regenerated and checked on restart
under an exclusive lock held through receipt creation.

The entire 2,304-pair/48-group/24-sequence-group software grid passed using
actual raw-MSA/PDB/SVD synthetic sources; all 11 false pair exports and six false
group exports were rejected, and recovery checked two committed triples.
Production source-closure I/O is stubbed only in software fixtures. Resources
were estimated before launch: two CPU/32 GiB/no swap/one BLAS thread, 32 GiB
output allowance, 100 GiB reserve, no GPU or new charges. Context/parent linkage,
shared mask/method/core qualification, domain/PAE/predictor controls and calibrated
phylogenetic/duplication inference remain required. Numerical signs do not infer
ancestral polarity or biological significance. [Full comparison workflow](full-triad-sequence-correspondence-20261002.md#full-correspondence-comparison-queued-october-2-utc).

### Full correspondence linkage to original contexts (October 2 UTC)

Full context linkage is implemented and queued after complete correspondence
comparison closure. Every original 283,409 context/566,818 designs/214,461 tied
references/428,922 sides retains its guide, family/taxon/gene/node/source/model
identity and lexical/tie order. Original parent, distinct model/version, three
pair-design and native-own/both-guide readiness are retained; shared physical
results do not promote excluded logical contexts. Missing results remain NULL,
separate from measured zero-pass bitmaps. Reference links retain eight normalized
physical keys for complete numerical/status/overlap access.

All 27 logical mask/method/core combinations include joint both-mask, both-method
and both-core requirements. The corresponding 48-order-pair bitmaps intersect;
full qualification requires all bits in each included leaf. Five lexical/any/all
tied-reference policies retain every missing/excluded/tied state; empty tie sets
fail the all-tie policy. Full scope includes 34,742,682 reference-screen cells,
91,824,516 context/design/scenario/screen cells, 459,122,580 repeated policy
decisions and 3,240 guide/design/scenario/screen/policy summary rows. These are
dependent eligibility screens, not independent events or statistical tests.
Conditional source/geometry qualification does not establish contrast-direction
stability, accepted duplication/phylogeny or biological significance.

The independent reader rebuilds ordered model SHA keys and every source gate
from raw context records, derives scenario bitmaps via unions of missing bits,
and reconstructs all policies, counts and denominators using SQL. Every source
record, output field/type, complete row grid and deterministic checkpoint must
agree before full source/hash and both original-journal closure. Interrupted
restart regenerates and checks 1,000-context chunks under an exclusive lock.
The full 28-context/60-tie/27-scenario software grid passed 7,560 policy decisions,
540 one-screen summaries and recovery of a four-context chunk; all 18 false
exports were rejected. Production source-closure I/O is stubbed only in fixtures.
Resources were estimated before launch: two CPU/32 GiB/no swap/one BLAS thread,
32 GiB output allowance, 100 GiB reserve, no GPU or new charges. Production and
prediction/domain/PAE/phylogenetic dependence/sampling/calibration remain pending.
[Full context workflow](full-triad-sequence-correspondence-20261002.md#original-context-linkage-and-joint-qualification-queued-october-2-utc).

## Full structural contrast sensitivity added October 2 UTC

The full closed set of 27,056 source-ready ordered physical triples is analyzed
under nine full/pLDDT70/both-mask by reference-common/cycle-consistent/both-core
scenarios, retaining all eight orders in each selected leaf. This creates
243,504 contrast groups and 1,461,024 decisions under six existing quality
screens. Signed RMSD(A, reference) minus RMSD(B, reference) extrema retain
original model versions and A/B/reference roles. Available/expected counts and
unique-fit status distinguish missing, nonunique, positive, negative, near-zero
and sign-uncertain dispositions. Both exact-zero and 1e-9 Å numerical boundaries
are exported; these thresholds are numerical classification boundaries and
do not define biological effect size. Sensitivity envelopes are not uncertainty
intervals. Eligibility requires every selected order to pass combined core and
inherited three-pair screens with original protein-length coverage denominators.
Source logical parents, native guides, reference choices/ties and family/taxon
dependence remain necessary gates for subsequent comparative inference.

The producer combines closed all-order summaries; a separate raw-fit SQLite
reader reconstructs every group, source role/grid, bitmap, status, range, sign,
count, denominator and deterministic checkpoint from all 865,792 original fit
rows. Full hash verification and both original completion/resource journals gate
closure. Software contracts passed 384 raw rows/108 scenario groups/216 decisions,
committed checkpoint recovery and 16 false-export rejections. Full production
and independent reconstruction passed all 243,504 groups/1,461,024 decisions;
two-original-journal closure binds 464,919 hashes. At 50 residues/70% coverage,
the joint-mask/joint-core quality cohort contains 5,448 physical triples:
1,975 uniformly positive, 1,785 uniformly negative and 1,688 sign uncertain
under the numerical boundary. These are dependent descriptive measurements;
context/tie gates and biological/statistical calibration remain required.
[Scope, results and resources](full-triad-contrast-sensitivity-20261002.md).

## Full original-context contrast directions added October 2 UTC

The closed physical contrast-sensitivity results are linked to all 283,409
unchanged original gene-tree contexts, retaining 214,461 tied references and
428,922 logical sides. Both native guides and availability/sequence-first
designs preserve original family/taxon/gene/node/model/version/source identities,
reference order and parent/distinct-model/three-pair/native gates. All nine
mask/core scenarios and six screens create 11,580,894 reference cells,
30,608,172 context/design cells and 153,040,860 repeated five-policy decisions.
These are dependent descriptions, not independent biological events or tests.

Lexical references remain fixed under each of three source gates. Any-eligible
native-both policy direction uses every quality-eligible tie and retains
eligible counts and four-category direction support; it does not stand in for
complete original-tie agreement. All-tie policies require every original tie
eligible and a nonempty set. Differing reference direction categories remain
disagreement, missing/near-zero/excluded states remain distinct, and geometric
eligibility does not establish directional consistency. Physical extrema and
1e-9 Å numerical direction classes remain descriptive, not confidence intervals
or biological effect cutoffs. Independent raw-model SHA/source-gate reconstruction
and SQL lexical/pool/support/count/denominator checks cover all original records
and deterministic checkpoints, followed by full hashes and two original journals.
Software checks passed 24 contexts/52 ties/2,160 decisions, committed checkpoint
recovery and 20 false-export rejections. Full production and independent readback
passed all 153,040,860 policy decisions and 10,800 count rows; final closure
verifies 465,222 bindings and both original completion/resource journals.
At ≥50 residues/≥70% original coverage, both masks/cores and native-both fixed
lexical references qualify 5,365 contexts per availability guide and 3,764/3,765
sequence-first contexts. Their respective positive/negative/sign-uncertain
counts are profile availability 1,951/1,755/1,659, MAFFT availability
1,950/1,754/1,661, profile sequence first 1,401/1,238/1,125, and MAFFT sequence
first 1,400/1,238/1,127. Complete any/all-tie counts and full exclusions are
published separately. Overlapping conditional cohorts cannot be pooled as
independent events. Figures passed exact-value checks and actual PNG/rendered
PDF review after resolving annotation overlap. Phylogenetic/predictor/domain/PAE
and calibrated inference remain required.
[Complete context semantics and resources](full-triad-context-contrasts-20261002.md).

## Joint sequence/structural contrast directions queued October 2 UTC

The complete physical correspondence control combines all 649,344
sequence-derived and 865,792 structure-derived raw fit dispositions across
27,056 original source-ready ordered triples. All 27 mask/method/core scenarios
retain separate source and joint signed AR-minus-BR RMSD extrema, expected and
available orders, uniqueness, strict-zero and 1e-9 Å numerical direction states.
Selected source geometry keys are deduplicated: the most inclusive joint
envelope contains 24 sequence plus 32 structural fits, not 384 repeated paired
comparisons. Matching source direction categories can include matching sign
uncertainty and do not establish stable sign or effect concordance.

Each screen intersects all selected original 48-bit six-by-eight order-pair
maps across masks, methods and cores; every bit must pass with inherited
three-pair exclusions and original protein-length denominators. No favorable
method/core/order or mask/method diagonal can create qualification. Full scope
is 730,512 joint groups, 4,383,072 dependent six-screen decisions and 1,134
disjoint qualified-direction count rows. Raw-fit SQLite reconstruction
independently checks complete grids/model roles/statuses/extrema and separately
derives Cartesian order-pair flags and complement-union qualification; all
export fields/types, counts, denominators and checkpoints require agreement.
Full hashes and both original journals gate completion.

Software contracts passed 14 triples/784 raw fits/378 joint groups/756 decisions,
committed two-triple recovery and completed-restart refusal; all 20 rehashed
false exports were rejected. Production source-closure I/O is stubbed only in
fixtures. Two CPU/32 GiB/no swap/one BLAS thread, 16 GiB output and 100 GiB
reserve were estimated before the three original jobs were queued. Production
awaits full correspondence-comparison closure. Sequence mappings use the same
predicted coordinates; predictor circularity, original contextual parent/native/
reference/tie gates, matched backgrounds, phylogeny and calibrated inference
remain required. All eight aims remain incomplete; GPU prediction stays paused.
[Complete joint-direction workflow](full-triad-joint-directions-20261002.md).


## Full contextual joint directions queued October 2 UTC

All 283,409 unchanged contexts/566,818 designs/214,461 tied references/428,922
sides are linked to the complete joint physical direction controls under all
27 mask/method/core scenarios and six screens. This defines 34,742,682 reference
cells, 91,824,516 context/design cells, 459,122,580 repeated five-policy direction
decisions and 32,400 disjoint state-count rows. Original guide/family/taxon/gene/
node/model/version/reference-order identities and parent/distinct-model/three-pair/
native gates remain unchanged. Normalized keys retain explicit mask, method and
core axes. Missing/gated direction remains null, distinct from measured near-zero.

Fixed lexical policies retain each of three source gates. Any-eligible native-
both direction uses every eligible tie with gate counts and four-category
support; agreement among incomplete subsets does not establish complete
original-tie agreement. All-original policies require nonempty complete
eligibility. Reference-category disagreement remains separate from geometry
quality; favorable reference/method/core/order choices cannot promote a sign.
The inherited numerical threshold does not define effect size or uncertainty.

Independent original-role SHA and raw context-gate reconstruction plus SQLite
lexical/all-pool/support/denominator checks require exact export fields/types,
source row order and all deterministic checkpoints, followed by full hashes
and both original invocation-linked completion/resource journals. Software
contracts passed 30 contexts/64 ties/216 joint physical groups/all 27 scenarios/
8,100 decisions, committed four-context recovery and completed-restart refusal;
all 23 false exports were rejected. Production source-closure I/O is stubbed
only in fixtures. The initial fixture receipt write failed on a missing output
directory; directory creation was fixed and the complete suite rerun to a new
v2 identity before launch. Two CPU/32 GiB/no swap/one BLAS thread, 192 GiB output
and 100 GiB reserve were estimated before three original queued launches.
Full production awaits physical joint-direction closure. Predicted coordinates
remain sequence-derived, contexts/settings remain dependent, and matched
backgrounds/predictor/domain/PAE/accepted phylogeny/reconciliation/calibration
remain required. All eight aims incomplete; GPU prediction paused.
[Full integration, semantics and reproduction](full-triad-context-joint-directions-20261002.md).


## Full whole-protein optimization provenance handoff queued October 2 UTC

The final handoff links the entire 75,070-input/five-tree/375,350-fit grid to all
75,070 original five-fit audit shards and every flagged follow-up. Each original
input/fit checksum, audit identity/five-fit checksum set/result digest and exact
manifest result is required. Follow-up scope equals the entire original census;
producer/reader keys equal all original flags. Original and follow-up errors
retain unverified numerical status, and reviewed outcomes cannot become passed
fits through hashing alone. Four original fit/audit/follow-up/reader processes
must be absent with invocation-linked original completion/resource journals and
full hashes before a final completion locator can be written. Candidate numerical
likelihood/coefficient/covariance/gradient replay remains the existing full
reader's responsibility; this handoff binds its exhaustive evidence. Software
contracts passed 15 fixture original audit states, all pass/review/error paths
and 11 rehashed corrupted-evidence rejections. Prelaunch estimates specify two
CPU/32 GiB/no swap, 4 GiB archives and 100 GiB reserve. Full native results,
versioned complete comparison integration and calibrated inference remain
pending; all eight aims incomplete, GPU prediction paused.
[Full scope, evidence and reproduction](full-whole-protein-optimization-closure-20261002.md).

## Full whole-protein comparison integration queued October 2 UTC

The new immutable export integrates every verified dynamic optimization
follow-up over all 75,070 frozen numerical inputs/five tree alternatives/
375,350 fits, preserving all 4,147,200 original comparison links. Original
passes/errors remain unchanged. Passing closed follow-ups select their fitted
refinement/recovery; unresolved follow-ups retain explicit review. Follow-up
errors retain the original audited flagged candidate with both source identities.
Selected coefficients, covariate scales, scaled/raw conditional GLS covariance,
variance ratios and residual scale are exported in a separate normalized
375,350-record parameter file. Conditional covariance is not a calibrated
interval or propagation of phylogenetic/predictor/sampling/model uncertainty.

An independent selection implementation reconstructs every source summary and
parameter, validates dimensions, finite values and unit transformations, then
uses the existing separate likelihood arithmetic/decision reader for every link.
Different observations, original errors, optimization/support reviews and
negative nested gains remain explicit. Deterministic parameter compression and
source-bound tree checkpoints permit identical-plan recovery; completed exports
refuse restart. Software checks passed all five relationships and selection
paths, committed recovery, completed-restart refusal and 26 false-source/export
rejections. Full production requires prior four-journal optimization closure,
every source/artifact hash and final two-original-journal closure. Prelaunch
resources: two CPU/32 GiB/no swap/one BLAS thread, 64 GiB output and 100 GiB
reserve; no GPU/new charges. Full biological/model calibration remains required.
The older five-candidate comparison version remains separate, and these older
measurement fits do not establish expanded-atlas coverage. All eight aims remain
incomplete; GPU prediction remains paused.
[Full source-selection semantics and reproduction](full-whole-protein-comparisons-20261002.md).

## Expanded background and correspondence controls completed October 2 UTC

Full expanded background scope closed at 146,172 physical pairs/584,688
directed states, preserving old/new source ownership and all original missing,
short, nonunique and RMSD flags. Both-order original-length screens passed
complete independent reconstruction. At 50 residues/70% original coverage,
59,844 background pairs pass both masks; this is not a matched effect estimate.
All 4,250,692 fixed metadata selections were carried through coverage without
rematching, including unmatched and excluded records. Independent numerical
balance/reuse reconstruction closed 62,208 feature rows/7,776 coverage strata/
32,682,096 reuse rows. Across all 54 fixed scenarios, guide/policy/mask/screen/
feature groups retain signed extrema, absolute medians/maxima and unavailable
counts. The complete 4,608-row summary was independently recomputed. Figures
use all both-mask/n50c70 group maxima, retain distinct SD denominators and do
not establish universal balance thresholds. Maximum matched target–control
SMD0.1676 coexists with retained-target shifts up to1.2182 baseline SD relative
to all targets. Matched balance does not qualify representative event coverage.
[Full results and inspected figures](full-screened-background-results-20261002.md).

All 324,672 native sequence alignment states, 649,344 sequence geometry
dispositions and 10,389,504 sequence/structural order-pair comparisons also
passed independent readback and source/journal closure. The complete joint
physical-direction control retains all 1,515,136 raw fit states and 730,512
mask/method/core groups. In the both-mask/both-method/both-core n50c70 cohort,
5,402 triples pass quality:1,416 positive/1,270 negative/2,716 sign uncertain.
Numerical direction envelopes are sensitivity descriptions, not calibrated
uncertainty or evolutionary polarity. Original-context eligibility linkage is
closed; contextual joint-direction readback/closure remains pending. Both
figures passed full count/SVG/PDF-label checks and actual PNG/rendered-PDF
inspection. Complete new source hashes were refreshed alongside original
terminal journals in runtimev20 (2,301,852 distinct bindings). Expanded
working-model inputs, predictor/domain/PAE, accepted phylogeny/reconciliation,
sampling/dependence and calibrated effects remain required. All eight aims
incomplete; GPU prediction paused.


## Full first-horizon ancestral accounting, October 2

All 1,620 initial independent-chain dispositions were source-closed against
original configurations, input hashes, native output/sample-audit/scalar-log
artifacts, diagnostic reports and both original completion/resource journals.
The 1,000-iteration horizon produced 1,617 integrity-checked chains and three
native memory-allocation failures. Diagnostics covered 402 four-chain quartets
at burn-in cutoffs 250 and 500, with no quartet passing every scalar check.
Allocation notices were censused from all original stderr files and retained
for five chains, including three integrity-passing outputs. These are source
and diagnostic accounting results, not joint posterior qualification.

Three failed chains receive separate same-seed attempts with only the native
address-space limit changed from 12 to 48 GiB. Serial execution uses one CPU,
64 GiB total memory and no swap. Original attempts remain unchanged, with no
sample concatenation. Corrected quartet lineage, new integrity/diagnostic
checks and length/alignment/categorical/scalar mixing qualification remain
required. See [complete sources and limitations](baliphy-method-assessment-20260927.md#full-initial-horizon-completed-and-source-closed-october-2).


## Full ancestral diagnostics and initialization census, October 2

State extraction and length/categorical diagnostic accounting closed all
1,617 processed chains, 402 complete quartets and three unresolved native
failures, with 49,902 source/artifact hashes and three original controller
journals. The complete 810-row publication retains every group/cutoff and
all count partitions; no quartet passes every scalar or candidate-length
screen. Constant/no-observed variation states retain review status. The
8-panel PNG/SVG/PDF figure passed source-count, table, SVG-metadata, PDF-value
and visual checks. Anchors and temporal patterns remain correlated units.

A separate full original native-log census checked 69,601,090 numeric cells
across all 1,620 chains and exported nine declared scalar ranges per chain.
All 335 nonfinite gamma-shape cells occurred in iterations 1–14, before both
original burn-in cutoffs. Matching installed v4.3 source explicitly handles
infinite gamma shape using unit rates; displayed infinity alone is not
treated as proof of a defective model. The first higher-memory isolated
recovery still failed, so resource allowance is not treated as sufficient
remediation. Priors, taxa, sequences and scientific thresholds are unchanged.
[Full sources, figures, counts and interpretation](baliphy-full-first-horizon-diagnostics-20261002.md).
