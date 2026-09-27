# Independent-control coverage for whole/domain cases

None of the 39 exact case sequences occurs in the frozen 25,322-model production
ESMFold inventory or the 79 ESMFold experimental-control models across four
length tiers. None matches either target or experimental-entity sequence hashes
in the existing 3,592-row experimental sequence-screen inventory. These are
negative results for the declared local inventories, not evidence that no
experimental homolog exists.

A new experimental-filtered RCSB accession query covered all 39 unique case
accessions and returned zero polymer entities. The exact accession list, query,
raw response and retrieval timestamp are cached. Readback checked query membership,
experimental-only filtering, full result counts and source/response hashes. The
query nominates entities by database accession; related sequences, constructs or
unlinked structures remain outside this search.

`audit_case_independent_prediction_coverage.py` prepared a FASTA containing all
39 sequences (19,561 residues), reconstructed from the verified full-coordinate
inputs and checked against the case sequence hashes. This is a future input,
not a prediction job. No GPU work was started. The project-level GPU pause still
applies; the new case set is distinct from the previously reported missing marker
sequences and must not silently be added to that count.

The immediate scientific consequence is that the observed coordinate contrasts
and their PAE qualification lack an existing exact-sequence independent control
in these inventories. Next options are a documented sequence-homology search for
experimental comparators, with coverage and construct checks, or independent
predictions when GPU use is resumed. Homologous experimental proteins would be
contextual controls, not exact validation of these fungal duplicate structures.

## Reproduction

1. Run `scripts/audit_case_independent_prediction_coverage.py` to check all frozen
   inventory hashes and prepare the exact-sequence snapshot.
2. Run `scripts/inventory_experimental_accession_matches.py --mapping
   results/structural_comparisons/case-independent-control-inputs-20260927-v1
   --output <new-output-directory>` to repeat accession nomination.

Local tables and FASTA are in
`results/structural_comparisons/case-independent-control-inputs-20260927-v1`;
RCSB results are in
`results/experimental_structures/whole-domain-case-accessions-20260927-v1`.
Versioned provenance and readback are in
`metadata/case_independent_control_coverage_completed_20260927.json`.
The accession query used one HTTP worker, with no paid resource or GPU use.

## Experimental sequence-similarity search

The subsequent search completed all 39 full sequences using RCSB's protein
sequence service, an identity cutoff of 0.30 and E-value cutoff of 1e-5, restricted
to experimental entities. Every hit was retained, without representative
selection or a coverage filter. The method follows the
[RCSB Search API](https://search.rcsb.org/) query schema; these thresholds nominate
comparators and do not establish functional equivalence.

Twenty-seven sequences had hits and twelve had none at these settings. The
results comprise 3,091 sequence/entity alignments and 1,205 unique experimental
polymer entities. Raw verbose responses preserve reported alignments and scores.
Every query and response count passed readback. Each query alignment was
reconstructed against its exact original sequence; aligned-column and paired-
residue counts, identity denominators, and query/subject span coverage are
exported explicitly. Subject spans have valid lengths. Subsequent independent subject-sequence
verification passed for every alignment; coordinate verification remains separate. Query alignment
coverage must not be reported as experimentally observed residue coverage.

Run `scripts/search_case_experimental_homologs.py` for the fixed restartable
search and `scripts/readback_case_experimental_search.py` for the full readback.
The former caches each request/response under a fixed configuration and preserves
all 39 query dispositions, including no-hit queries. The latter waits on the
recorded process identity and requires successful termination. Both services
completed successfully without GPU or paid-resource use.

Raw search results are in
`results/experimental_structures/whole-domain-case-sequence-search-20260927-v1`;
validated alignment-context rows are in
`results/experimental_structures/whole-domain-case-sequence-readback-20260927-v1`.
Versioned hashes and counts are in
`metadata/case_experimental_sequence_search_completed_20260927.json`.
Next, entity sequences, constructs, domain coverage, experimental quality,
observed residues and model-training overlap must be checked before selecting
structural controls. Results from related proteins are not exact-sequence
validation of the fungal cases.

## Complete entity and entry metadata verified

All 1,912 metadata requests completed for the 1,205 entities and 707 entries.
Both retrieval and independent subject readback terminated successfully. All
3,091 search alignment subjects exactly match the retrieved canonical entity
sequences. Response identities, hashes, mutation annotations, experimental
methods, resolution and release dates are retained. The completion record is
`metadata/case_experimental_subjects_completed_20260927.json`.

The metadata classifies 706 entries as experimental and one as integrative.
`retrieve_case_experimental_coordinates.py` downloads all 706 experimental
entries, retaining the integrative entry explicitly as an exclusion. Selection
does not depend on quality or agreement with predictions. Retrieval uses one
HTTP worker, one CPU, 2 GiB memory and no swap; disk allowance is 10 GiB with
20 GiB minimum headroom. The planning allowance was 0.5–6 hours. No GPU or
paid infrastructure is used.

`scripts/readback_case_experimental_archives.py` waits for the recorded producer
process and successful terminal state, then independently streams every archive
to verify identity, hashes, compressed/uncompressed sizes and gzip integrity.
It reconstructs methodology selection from the metadata. The queued
`scripts/map_case_experimental_ca_residues.py` requires the successful archive
readback before mapping all nominated entities, chains, deposited models and
full canonical sequence positions to CA observations. Missing residues,
alternates, partial occupancy and modified/mismatching monomers remain explicit.
It checks coordinate entity sequences against the retrieved experimental
sequences, not against the fungal targets. It uses one CPU, 4 GiB memory,
no swap and a 10 GiB output allowance, with 1–12 hours planned.

Metadata and subject-readback directories are respectively
`results/experimental_structures/whole-domain-case-metadata-20260927-v1` and
`results/experimental_structures/whole-domain-case-subject-readback-20260927-v1`.
Coordinates, archive checks and residue maps are under the same parent with
suffixes `whole-domain-case-coordinates-20260927-v1`,
`whole-domain-case-coordinate-readback-20260927-v1` and
`whole-domain-case-ca-mapping-20260927-v1`. Launch records preserve commands,
process identities, script hashes and resource limits in `metadata/`.

Canonical sequence agreement does not establish observed coordinate coverage,
construct equivalence, training independence or experimental quality. Mapping
outputs will need full readback and projection through the query/subject
alignments before applying domain/outside coverage thresholds.

## Domain and outside-region alignment coverage

`scripts/screen_case_experimental_domain_coverage.py` reconstructed paired query
residue positions in every one of the 3,091 search alignments. Query residues
aligned to a subject gap do not contribute. Each case/domain boundary is screened
against its original interval length and, separately, original protein length
minus that interval. A shared experimental entity must pass the screen for every
A/B/reference interval and boundary. Separate alignment contexts are not combined
to manufacture residue coverage. All six screens and all thirteen cases, including
zero-candidate cases, are retained.

At n50/c70, shared entity counts are:

| Case | Domain coverage | Domain and outside coverage |
| --- | ---: | ---: |
| Heliocybe OG0000054 | 31 | 19 |
| Rhodonia OG0000095 | 0 | 0 |
| Cryoendolithus OG0000107 | 97 | 97 |
| Scedosporium OG0000152 | 23 | 0 |
| Furculomyces OG0000230 | 10 | 6 |
| Jaapia OG0000294 | 102 | 102 |
| Synchytrium OG0000336 | 0 | 0 |
| Neolecta OG0000972 | 3 | 0 |
| Leucosporidium OG0001082 | 0 | 0 |
| Smittium OG0001200 | 0 | 0 |
| Leucosporidium OG0001203 | 0 | 0 |
| Piloderma OG0002650 | 0 | 0 |
| Phycomyces OG0002812 | 26 | 26 |

These are candidate entity counts, not independent experiments or observed
coordinate coverage. Multiple entities may belong to the same PDB entry; entries
may share constructs or experimental context. Both the fungal prediction coverage
and experimental alignment coverage are needed: Cryoendolithus and Furculomyces
passing this search screen does not override their earlier predicted-coordinate
coverage limitations. The metadata retrieval and subject-sequence check have completed; coordinate
and observed-residue validation are now underway.

Full alignment/interval decisions, shared-entity decisions and all 78 case/screen
rows are in
`results/experimental_structures/whole-domain-case-domain-coverage-20260927-v1`.
Serialized readback and source hashes are recorded in
`metadata/case_experimental_domain_coverage_completed_20260927.json`.

## Archive retrieval and full readback completed

All 706 experimental archives downloaded successfully (962,111,990 compressed
bytes; 4,392,908,993 uncompressed bytes). Retrieval and independent archive
readback both terminated successfully. Every archive passed identity, hash,
size and complete decompression/CRC checks; all plaintext hashes are retained.
The closure is `metadata/case_experimental_coordinates_completed_20260927.json`.
The residue-mapping service has now passed its prerequisite gate and is writing
per-entry residue tables. Those tables are not yet a completed or independently
validated coordinate-coverage analysis.

## Full raw-coordinate readback queued

`readback_case_experimental_ca_mapping.py` waits on the exact recorded mapping
PID, creation time and command, then requires successful terminal service state.
It checks all 706 entries against the original mmCIF files, not a sample: every
selected entity sequence, chain/model grid, full sequence position, exported CA
atom field and observation status. It also checks that no selected raw CA record
was omitted. Status logic is separately implemented; eight fixtures covering
missing, multiple, modified, invalid, alternate/partial and unambiguous records
passed. The producer and checker share Biopython's mmCIF parser, so parser
independence is not claimed.

The checker will export `chain_model_coverage.tsv` under
`results/experimental_structures/whole-domain-case-ca-readback-20260927-v1`.
This is experimental entity coverage, not fungal-query or domain coverage;
projection through all recorded search alignments remains required. The launch
record is `metadata/case_experimental_ca_readback_launch_20260927.json`.
Resources are one CPU, 4 GiB RAM, no swap, 0.2 GiB output and a 1–12 hour planning
allowance. The checker is queued; its validation has not yet completed.

## Observed domain/outside coverage queued

`scripts/screen_case_experimental_observed_coverage.py` is queued behind the
exact full-CA-readback process and successful terminal state. It reconstructs
both query and subject positions in every search alignment, then projects only
unambiguous, standard-monomer, full-occupancy CA positions onto the fungal query.
Alignment offsets, insertion/deletion gaps, observed-position projection and
invalid double-gap rejection passed fixtures recorded in
`metadata/case_experimental_observed_coverage_fixture_checks_20260927.json`.

Every case, role, interval boundary and all six coverage screens are retained.
A single experimental chain/model must qualify for every required role/interval;
residues from different chains, models or alignment contexts are not combined.
The integrative entry is excluded explicitly. Denominators remain the original
fungal domain and outside lengths, not the number of resolved experimental
residues. Summary counts distinguish chain/models, entities and PDB entries;
none is asserted to represent independent experiments.

The output is planned under
`results/experimental_structures/whole-domain-case-observed-coverage-20260927-v1`.
One CPU, 8 GiB RAM, no swap and 2 GiB output are allocated, with a 0.5–6 hour
planning allowance. `metadata/case_experimental_observed_coverage_launch_20260927.json`
records the exact process and script hash. The service is queued, not a completed
coverage result. Geometry, construct quality and prediction-training overlap
remain separate checks.

## Method, construct and starting-model annotations reviewed

`summarize_case_experimental_metadata.py` completed an annotation review of all
707 entries and 1,205 candidate protein entities. It retains method-specific
refinement/validation fields, resolution arrays, release dates, primary citations,
source/host organisms, construct annotations and all reported starting models.
No candidates were ranked or filtered on these fields. All serialized outputs
passed readback; all release dates and entity mutation/nonstandard counts and
mutation text were checked again against raw metadata.

Six entries explicitly mention AlphaFold starting models: 8S0B, 8S0D, 8S0E,
9TI8, 9TI9 and 9U4X. Six others report other computational starting models;
405 report only experimental starting models, 283 lack this annotation and
seven have other/incomplete annotations. These categories cover the full
707-entry inventory, including its one integrative entry. Target-chain
attribution remains unresolved. Missing annotation does not establish absence
of prediction-assisted refinement, and experimental starting models may have
their own dependencies.

There are 117 entities with positive reported mutation counts. Zero counts are
not proof that constructs match native fungal proteins. All candidate homologs
still require sequence, coordinate and biological-context review. Release dates
are recorded but are not treated as prediction training/template cutoffs;
training independence remains unresolved for every entry.

Outputs are `entries.tsv`, `entities.tsv` and
`method_specific_annotations.json` in
`results/experimental_structures/whole-domain-case-metadata-review-20260927-v1`.
Hashes, counts and explicit flagged entry IDs are recorded in
`metadata/case_experimental_metadata_review_completed_20260927.json`.

## Case-specific metadata linkage completed

`link_case_experimental_metadata.py` joined all 7,230 candidate entity/screen
rows to the exact entity and entry annotations. It retained all 156 combinations
of thirteen cases, six screens and two coverage regions, including zeros. An
independent pandas join recomputed every summary count and flagged-entry list
from the original tables; all matched.

At n50/c70 for **alignment** coverage of both domain and outside regions:

| Case | Experimental entities | PDB entries | Entries with AlphaFold starting annotation | Entities reporting mutations |
| --- | ---: | ---: | ---: | ---: |
| Heliocybe OG0000054 | 19 | 19 | 0 | 1 |
| Cryoendolithus OG0000107 | 97 | 97 | 3 | 0 |
| Furculomyces OG0000230 | 6 | 6 | 0 | 2 |
| Jaapia OG0000294 | 101 | 100 | 0 | 3 |
| Phycomyces OG0002812 | 26 | 26 | 0 | 17 |

The other eight cases have zero candidates for both regions at this setting.
The earlier Jaapia count of 102 includes the integrative entry; methodology
stratification now distinguishes 101 experimental entities and one integrative
entity. No candidate was discarded because of mutation or starting-model flags.
These annotations require target-chain and construct review, and absence of a
flag is not independence. The table does not establish resolved-coordinate
coverage, nor does it override fungal-prediction coverage requirements.

Artifacts are `annotated_candidate_entities.tsv` and
`case_metadata_summary.tsv` under
`results/experimental_structures/whole-domain-case-metadata-links-20260927-v1`.
Full source bindings and the independent readback are recorded in
`metadata/case_experimental_metadata_links_completed_20260927.json`.

## Exact query/experimental residue correspondences completed

`prepare_case_experimental_residue_pairs.py` froze all 3,091 alignment contexts
for the 39 queries. It records 1,419,058 paired residue occurrences, including
621,030 identical pairs, plus 138,892 query residues and 68,232 subject residues
aligned to gaps. All positions are one-based in the original canonical sequence.
Each mapping was reconstructed against both exact sequences and independently
checked using alignment-column indexing. All twelve no-hit queries remain in
the sequence-disposition table.

These mappings are inputs for intersecting fungal A/B/reference correspondences
with a single observed experimental chain/model. They do not themselves establish
observed coverage, equivalent conformational states or a common-residue geometric
fit. The integrative candidate remains in this sequence-only inventory and must
remain excluded from experimental-coordinate analyses. All occurrences and
contexts are retained rather than counted as independent observations.

Output: `results/experimental_structures/whole-domain-case-residue-pairs-20260927-v1`.
Closure: `metadata/case_experimental_residue_pairs_completed_20260927.json`.

## Shared four-protein sequence correspondences completed

`prepare_case_experimental_quartets.py` intersected the three fungal-to-
experimental residue maps for each of the thirteen distinct A/B/reference
model triplets. Eight triplets have at least one experimental entity shared by
all three sequences. All 923 shared-hit/context combinations are retained,
comprising 433,680 common residue occurrences; five triplets have explicit
no-shared-candidate dispositions. Both species-guide links remain attached.

Each row lists one-based positions in A, B, the fungal reference and the same
experimental entity. Every intersection passed an independent ordered membership
scan. All alignment-context combinations are retained without merging their
residues. These maps define sequence-anchored correspondences, distinct from the
structural-alignment maps used in earlier within-fungal comparisons. This
alternative definition must remain explicit in downstream sensitivity analyses.

No coverage or geometric acceptance follows from a shared sequence hit. Before
fitting, all four coordinates must be present, the experimental entry must be
classified experimental, the selected chain/model must have unambiguous CA
observations, and fungal confidence masks and domain/outside partitions must be
applied. Common-residue coverage uses original lengths, and degenerate fits must
remain excluded from rotation interpretation.

Outputs are in
`results/experimental_structures/whole-domain-case-sequence-quartets-20260927-v1`;
closure is `metadata/case_experimental_sequence_quartets_completed_20260927.json`.

## Observed quartet filtering and direct geometry queued

`fit_case_experimental_quartets.py` is queued behind the exact full experimental
CA readback and successful terminal state. It binds the 923 sequence quartet
maps to the exact fungal model versions and checked coordinate files. Every
experimental chain/model is retained, with full and pLDDT70 fungal masks.
Only positions observed in all four proteins contribute. The integrative entry
is explicitly excluded and all five no-shared-candidate triplets remain recorded.

For each of the 26 domain-boundary triplets, common quartets are partitioned into
inside all three fungal domains, outside all three, and mixed membership. Mixed
positions remain in the whole-protein comparison but are excluded from domain
and outside fits. All six protein pairs (AB, AR, BR, AE, BE, RE; E denotes the
experimental homolog) use identical four-way residue sets for each region.
Proper rotations are fitted separately to whole, domain and outside regions;
the domain transform is additionally evaluated on outside residues without
refitting. Insufficient or degenerate fits remain explicit. Original lengths
and common position sets are exported for subsequent coverage checks.

Four rigid-motion/outside-translation fixtures passed. Each computed production
fit will compare SVD and quaternion RMSDs, including outside residuals under the
domain transform. These inline checks do not replace a full independent readback.
Coverage qualification, construct/dependence review and interpretation remain
outstanding; these homolog comparisons are not exact fungal prediction controls.

The service uses one CPU, single-threaded BLAS, 8 GiB RAM and no swap, with a
5 GiB output and 1–12 hour planning allowance. No GPU use is involved. Output is
`results/experimental_structures/whole-domain-case-quartet-fits-20260927-v1`;
launch identity is `metadata/case_experimental_quartet_fits_launch_20260927.json`.
The stage is queued, not completed.

## Residue mapping production completed; full raw readback active

The mapping service terminated successfully after all 706 entries. Its output
contains 2,136 chain/model grids and 1,280,073 canonical-position rows. All
per-entry receipts and compressed table hashes passed verification, and all
reported category partitions sum to the aggregate totals:

| Production classification | Position rows |
| --- | ---: |
| Unambiguous full-occupancy CA | 1,027,254 |
| No CA observation | 244,982 |
| Multiple CA records | 3,072 |
| Alternate or partial occupancy | 3,188 |
| Nonstandard or mismatching monomer | 801 |
| Invalid coordinates or occupancy | 776 |

These are production classifications pending full independent raw-coordinate
readback. Counts include multiple chains and models and are not independent
residues or experiments. The readback has started processing the original files;
coverage and geometry services remain gated on its successful completion.
Production closure is
`metadata/case_experimental_ca_mapping_production_completed_20260927.json`.

## Experimental candidate repetition and publication links

`measure_case_experimental_dependence.py` checked all 1,204 experimentally
classified entities (the remaining entity is integrative). They contain 314
exact canonical sequences. A graph linking identical sequences, shared entries,
DOIs or PubMed identifiers has 196 components; the largest contains 541 entities.
All graph components matched an independently traversed adjacency graph. Missing
citation identifiers are not joined. Every candidate remains retained.

For n50/c70 alignment coverage of both regions:

| Case | Entities | Exact sequences | Linked components |
| --- | ---: | ---: | ---: |
| Heliocybe | 19 | 9 | 6 |
| Cryoendolithus | 97 | 5 | 4 |
| Furculomyces | 6 | 5 | 3 |
| Jaapia | 101 | 46 | 38 |
| Phycomyces | 26 | 10 | 6 |

All other cases retain zero counts at this setting. Components are formed on the
full candidate inventory before coverage selection, so links through other
entries remain represented. A publication or entry can connect different proteins;
these components are conservative dependence annotations, not proven independent
experiments, phylogenetic units or effective sample sizes. Neither 1,204 entities
nor 196 components should be treated automatically as independent replicates.
Geometry-based selection was not used.

Entity assignments, explicit linkage edges and all 156 case/screen/region rows
are in `results/experimental_structures/whole-domain-case-experimental-dependence-20260927-v1`.
The completion record is `metadata/case_experimental_dependence_completed_20260927.json`.

## Full observed-coverage readback queued

`readback_case_experimental_observed_coverage.py` waits on the exact coverage
producer and requires successful terminal state. It derives the complete expected
row set by expanding the prior alignment/interval/screen inventory over every
audited experimental chain/model. Every observed count is independently rebuilt
from the separately frozen residue-pair maps and audited CA positions, using
original domain/outside lengths. Exact rational arithmetic checks all thresholds.
It recomputes all role/interval qualifications and all 78 case/screen summaries,
including zero-candidate cases, and verifies distinct entity and entry counts.

The checker is queued with one CPU, 8 GiB RAM, no swap, 0.1 GiB output and a
0.25–4 hour planning allowance. Its output will be
`results/experimental_structures/whole-domain-case-observed-coverage-readback-20260927-v1`.
The launch identity is recorded in
`metadata/case_experimental_observed_readback_launch_20260927.json`.
This check is not yet completed and does not validate geometry or independence.

## Full quartet geometry readback queued

`readback_case_experimental_quartets.py` waits on the exact geometry producer and
successful terminal state. It reconstructs the complete chain/model/mask/domain/
context combination set, checks observed four-way positions against source
coordinate files, and rebuilds domain/outside/mixed partitions and denominators.
All eighteen fit rows per partition (six pairs × three regions) must be present.
Every RMSD and outside-under-domain-transform residual is recomputed using
separately written quaternion-matrix formulas; degenerate and insufficient fits
must have the correct status and blank metrics. The same numerical linear-algebra
library remains a shared dependency.

Twelve rigid/noisy fixtures and one collinear degeneracy fixture passed. Maximum
RMSD disagreement with direct SVD in those fixtures was 9.26e-16 Å. This is a
software check, not production validation. Production readback is queued with
one CPU, single-threaded BLAS, 8 GiB RAM, no swap, 0.1 GiB output and 1–12 hours
planned. Output will be
`results/experimental_structures/whole-domain-case-quartet-readback-20260927-v1`;
launch and fixture evidence are in
`metadata/case_experimental_quartet_readback_launch_20260927.json` and
`metadata/case_experimental_quartet_readback_fixtures_20260927.json`.
