# Methods draft: data preparation and exploratory structural comparisons

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
