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
Full and pLDDT70 C-alpha alignment inputs are being generated with original residue
positions retained. Input validation and expanded alignments remain pending;
this coordinate audit provides no PAE qualification or biological effect estimate.
