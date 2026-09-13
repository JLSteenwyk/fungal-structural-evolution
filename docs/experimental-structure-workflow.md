# Experimental structure benchmarking

Cross-predictor agreement cannot establish accuracy or remove shared
sequence-derived/training effects. Experimental references are therefore a
separate required benchmark, alongside predictor controls, confidence masking,
phylogenetic uncertainty and domain-orientation checks.

## Accession-linked candidate inventory

The first executable inventory queries all UniProt accessions associated with
the frozen 13,153 AlphaFold models, using the [RCSB Search API](https://search.rcsb.org/)
attribute `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession`.
It requests `polymer_entity` identifiers, all hits and experimental content only,
in deterministic batches of 100. This scope includes the full frozen marker-model
inventory rather than only the 266 predictor controls. Response bodies, queries,
HTTP statuses, retrieval times and checksums are preserved. HTTP 204 is recorded
as a completed no-result query; failures remain failures and are retried at most
four times. Result counts and uniqueness are checked before accepting a batch.

```bash
python scripts/inventory_experimental_accession_matches.py \
  --mapping results/structural_markers/gdm-expanded-v1 \
  --output data/experimental_structures/accession-inventory-v1
```

Existing batch responses are reusable only with matching configuration, queries
and checksums. A single worker uses a quarter-second pause after successful new
requests; the resource plan allows 2 GB RAM, 2 GB output and 0.1–3 hours on the
existing host. This stage does not download atomic structures.

Initial sequence-service connectivity checks returned a server error for an
obsolete parameter and timed out with the current documented sequence parameter.
The accession search and an atomic mmCIF endpoint responded successfully, so the
inventory proceeds through accession links. The documented sequence API remains
a possible later route for homologs lacking these links; no sequence-search
coverage is claimed from those connectivity checks.

## Candidate-to-benchmark gates

Accession association is candidate nomination only. Each returned entity needs
its deposited sequence, source organism, construct boundaries, mutations,
experimental method/quality and observed-residue mapping checked against our
exact protein. Missing/unobserved residues and engineered fragments must remain
explicit. Atomic coordinates need immutable source receipts and chain/entity
identity checks; experimental B factors are not predicted pLDDT scores.

Comparisons must distinguish monomeric predictions from complexes, bound ligands,
interfaces and alternative conformations. Deposit/release dates and predictor
training/template overlap require explicit review: agreement with a training-era
structure is not an independent generalization benchmark. Absence of an accession
hit does not establish absence of experimental homologs. None of these gates is
satisfied solely by completion of the inventory.


## Completed inventory

All 13,153 distinct model-associated accessions were queried in 132 batches,
yielding 3,592 distinct experimental polymer-entity candidates. Every saved batch
receipt/response hash, hit count and the complete deduplicated ID union passed
readback. Receipts and the candidate identifier list are versioned; full API
responses remain outside Git. This is a candidate count, not the number of
matched proteins, independent experiments or validated benchmark structures.
RCSB entity metadata was verified accessible for the next sequence/construct
screen; that screen has not yet been executed across these candidates.


## Candidate metadata retrieval executing

The complete inventory resolves to 3,592 polymer entities across 1,635 PDB
entries. `retrieve_experimental_metadata.py` retrieves full entity and parent
entry records through the [RCSB Data API](https://data.rcsb.org/), preserving
sequences, construct/organism annotations, accession mappings and entry-level
experimental/date metadata for subsequent screening. Two workers make 5,227
requests before retries, with a quarter-second pause after each new successful
request. The resource plan allows 2 GB RAM/output and 0.2–4 hours.

Each HTTP 200 response must parse and report the exact requested `rcsb_id`.
Atomic response/receipt writes record URL, retrieval time, configuration and
SHA256. Restarted retrieval checks cached bytes and identities; a whole-stage
completion receipt is written only after every requested record succeeds. Raw
metadata remain outside Git under `data/experimental_structures/metadata-v1`.
The initial entity responses were verified; the full retrieval is still running.

```bash
python scripts/retrieve_experimental_metadata.py \
  --inventory data/experimental_structures/accession-inventory-v1 \
  --output data/experimental_structures/metadata-v1
```

This is acquisition, not a sequence-matching or experimental-quality result.
Deposited polymer sequences and observed atomic residues will be assessed
separately; a complete deposited sequence may still have missing coordinates.


## Frozen initial deposited-sequence screen

The first explicit partial screen freezes 666 verified entity responses while
retrieval continues (2,926 candidates pending). Its 666 accession-linked model
comparisons comprise 541 exact full canonical sequences, 66 exact uniquely
placed fragments, 45 requiring alignment/variant review, 12 containing the target
inside a longer construct, and two with noncanonical sequences. These counts
follow retrieval order; they do not estimate benchmark coverage across all taxa.
Repeated experimental entities for one protein are not independent protein
observations. Exact-full rows passed sequence-hash and length equality readback.

`screen_experimental_sequences.py` retains unique fragment offsets and detects
repeated-fragment ambiguity without arbitrary alignment selection. It preserves
reported mutation, artifact and nonstandard-monomer counts, including unknowns.
Longer constructs are not automatically called affinity tags. Four tests cover
full matches, unique/repeated fragments, constructs, substitutions and missing
or noncanonical sequences. All benchmark-eligibility fields remain pending.

```bash
python scripts/screen_experimental_sequences.py \
  --metadata data/experimental_structures/metadata-v1 \
  --output results/experimental_structures/sequence-screen-partial-v1 \
  --allow-partial
```

The receipt pins every included entity response receipt. New runs against the
expanding download require new output directories and may include more entities.
Without `--allow-partial`, complete metadata retrieval is required. Even an exact
canonical sequence may include chemically modified monomers or missing atomic
coordinates; sequence equality alone does not establish benchmark eligibility.


## Atomic-coordinate acquisition for exact-sequence candidates

The frozen partial sequence screen nominates 215 distinct PDB entries containing
541 exact full canonical-sequence entity/model matches (34 distinct model
proteins). Every nominated entry is included in coordinate acquisition, without
selecting entries by observed agreement with predictions. `retrieve_experimental_coordinates.py`
streams the compressed mmCIF files from the RCSB coordinate service, checks full
gzip readability/CRC and the leading mmCIF data-block identity, then records
URL, time, sizes and SHA256. Two tests verify identity checking, complete-stream
reading and rejection of a truncated archive. These checks do not validate
atomic records or residue correspondence.

```bash
python scripts/retrieve_experimental_coordinates.py \
  --screen results/experimental_structures/sequence-screen-partial-v1 \
  --output data/experimental_structures/coordinates-exact-partial-v1
```

The single-worker resource plan allows 2 GB RAM, 20 GB disk headroom and 0.5–8
hours, retaining complete complexes. All 215 archives completed integrity and data-block checks. The expected entry
universe, per-entry receipts and every compressed-file hash passed readback. Subsequent screens must account for
multiple models, chains, alternate locations, occupancy, unresolved residues,
chemical modifications and complex context. Exact canonical sequence is not
permission to treat every atom or entity as an eligible benchmark. The full
metadata inventory continues independently; these downloads reflect a frozen
partial candidate set, not complete experimental coverage of the project.


## Observed experimental CA correspondence executing

`map_experimental_ca_residues.py` parses all 215 downloaded entries and checks
each nominated entity's coordinate-file canonical sequence against the exact
model sequence hash. CA mapping uses mmCIF entity, model, label-chain and
[label sequence IDs](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.label_seq_id.html),
not author numbering. Every full-sequence position is emitted for each selected
chain and deposited model. Missing observations, multiple CA records,
nonstandard/mismatching monomers, invalid coordinates/occupancy and alternate or
partial-occupancy observations remain distinct. No arbitrary alternate-location
choice is made. Raw atom records retain author identifiers, insertions and B
factors; B factors are not pLDDT.

Two tests cover absent, duplicated, alternate/partial, modified, nonfinite and
zero-occupancy observations. Initial entries are producing both unambiguous and
missing-coordinate rows; the complete mapping is still running. Per-entry output
hashes and configuration checks support restartable processing without changing
already completed results. The resource plan allows one CPU, 8 GB RAM, 5 GB
output and 0.2–4 hours for 1.30 GB uncompressed input.

```bash
python scripts/map_experimental_ca_residues.py \
  --screen results/experimental_structures/sequence-screen-partial-v1 \
  --coordinates data/experimental_structures/coordinates-exact-partial-v1 \
  --output results/experimental_structures/ca-mapping-partial-v1
```

Unambiguous full-occupancy CA correspondence is not a complete experimental
quality criterion. Non-CA atoms, experimental method/resolution, refinement
validation, biological assemblies, ligands, alternate conformations and
prediction-training overlap still require assessment. This calculation does
not yet compare experimental and predicted coordinates.


## Full metadata and sequence screen completed

All 5,227 entity/entry metadata responses completed. Independent readback verified
every response/receipt hash, declared identity and the exact expected entity and
parent-entry universes. The complete sequence screen covers all 3,592 candidates:
2,850 exact full canonical sequences, 187 exact unique fragments, 485 requiring
alignment/variant review, 52 longer constructs and 18 noncanonical sequences.
The full matches represent 80 distinct predicted proteins across 1,036 PDB
entries. These repeated-entity counts do not imply independent protein controls.
The earlier 215-entry coordinate snapshot remains a subset; expansion and
quality screening of the full set are still required.

## Initial coordinate mapping completed and audited

All 215 entries completed CA correspondence mapping: 157,670 model/chain/position
rows, including 129,095 unambiguous full-occupancy CA positions, 28,405 without
CA observations, 78 with multiple CA records, 25 with nonstandard/mismatching
monomers and 67 with alternate or partial occupancy. Counts are repeated
observations across structures/chains/models, not distinct sequence positions.

The audit verified every exported position grid, reconstructed sequence hash,
embedded atom identity and count partition. It independently read original
mmCIF atom fields for 2,823 positions across five deterministically selected
entries. Raw atom fields for other entries were not independently reparsed.
This establishes the stated correspondence checks, not an experimental-quality
pass or a completed structural-accuracy comparison.

```bash
python scripts/audit_experimental_ca_mapping.py \
  --mapping results/experimental_structures/ca-mapping-partial-v1 \
  --screen results/experimental_structures/sequence-screen-partial-v1 \
  --coordinates data/experimental_structures/coordinates-exact-partial-v1 \
  --output results/experimental_structures/ca-audit-partial-v1
```


## Full reference metadata and coordinate expansion

The frozen `reference-metadata-v2` inventory joins all 2,850 exact sequence rows
(80 predicted proteins) to 1,036 verified entry responses. Every entry's method,
methodology and resolution fields passed an independent source-field readback.
The initial v1 inventory is retained outside Git; v2 explicitly adds methodology.

The search's experimental-content filter also returned four entries whose entry
metadata classify them as **integrative**: 8ZZ2, 9A15, 9A16 and 9A17. These remain
in the accession/sequence inventory but are deferred from the experimental
coordinate benchmark. RCSB describes [8ZZ2](https://www.rcsb.org/structure/8ZZ2)
as an integrative structure; experimental restraints do not make all coordinates
independent experimental observations. This correction applies to the full
inventory interpretation, not to sequence correspondence.

The 1,032 remaining entries comprise 583 electron-microscopy, 448 X-ray and one
solution-NMR structure. Full source method-specific geometry, map, refinement,
starting-model and date fields are preserved in a hash-tracked JSON artifact.
Missing fields remain unknown. Entry metrics cannot establish target-chain
quality; release dates alone cannot prove training independence. All 1,036
entries are classified as heteromeric protein (477), protein/NA (555), or
homomeric protein (4), so complex and conformational context require review.

Coordinate expansion is running with one streaming worker, 2 GB memory planning
and 20 GB disk headroom, with a broad 0.5–8 hour forecast on the existing host.
It reuses 215 prior archives after checking collection/configuration/entry
receipts, checksums, gzip integrity and entry identity; 817 downloads remain in
the full requested universe. No new paid services are involved. The new output
is `data/experimental_structures/coordinates-exact-full-v1`.

Reproduce the metadata stage with:

```bash
python scripts/summarize_experimental_reference_metadata.py --screen results/experimental_structures/sequence-screen-full-v1 --metadata data/experimental_structures/metadata-v1 --output results/experimental_structures/reference-metadata-v2
python scripts/retrieve_experimental_coordinates.py --screen results/experimental_structures/sequence-screen-full-v1 --reference-metadata results/experimental_structures/reference-metadata-v2 --reuse data/experimental_structures/coordinates-exact-partial-v1 --output data/experimental_structures/coordinates-exact-full-v1
```

Use new output paths to regenerate immutable metadata. The downloader resumes
only with identical pinned configuration. Tests cover archive identity and
truncation, verified reuse without network access, integrative deferral, and
rejection of a tampered source archive. Three tests passed. Quality review,
expanded residue mapping and direct prediction comparisons remain pending.


## Initial direct prediction–experiment agreement

`compare_experimental_predictions.py` completed comparisons against all 215
entries in the original audited CA snapshot: 633 unique target/chain/deposited
model combinations, each evaluated at predicted focal pLDDT 0, 70 and 90.
Experimental CA must be unambiguous and full occupancy; at least 50 matched
positions and half the complete canonical sequence are required. All deposited
models and chains are retained, with no choice based on prediction agreement.
Experimental B factors are never interpreted as prediction confidence. This
stage does not apply a six-residue context or PAE filter.

There are 1,803 accepted and 96 excluded threshold rows. All 1,803 accepted
comparisons passed a second geometry calculation using SciPy rotation and
condensed distances, independent of the production geometry helper, at 1e-8
relative/absolute tolerance. Predicted CA sequence hashes and coordinate files,
experimental complete sequence grids and mapping receipts were verified.

| Predicted pLDDT | Accepted comparisons | Proteins | Entries | Comparison-weighted median RMSD (Å) | Median of protein medians (Å) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 626 | 33 | 214 | 0.900 | 3.235 |
| 70 | 617 | 32 | 214 | 0.874 | 1.860 |
| 90 | 560 | 18 | 206 | 0.631 | 0.766 |

The difference between weighting schemes is substantial. Per-protein summaries
limit the influence of proteins with many deposited structures; neither summary
makes the retrieval-order partial reference subset representative of fungi.
Within each protein, its chain/model observations still receive equal weight.

Matched-cohort summaries retain identical target/chain/model combinations at
baseline and the higher threshold. For the 617 comparisons accepted at pLDDT 70,
median RMSD is 0.892 Å before filtering and 0.874 Å after. For the 560 accepted at
90, the corresponding values are 0.827 and 0.631 Å. These compare different
residue subsets of the same chains and do not establish a causal improvement.
The full 1,899-row threshold grid and coverage eligibility passed readback.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/compare_experimental_predictions.py --mapping results/experimental_structures/ca-mapping-partial-v1 --screen results/experimental_structures/sequence-screen-partial-v1 --predictions results/structural_markers/gdm-expanded-v1 --references results/experimental_structures/reference-metadata-v2 --output results/experimental_structures/prediction-agreement-partial-v1
python scripts/summarize_experimental_agreement.py --comparisons results/experimental_structures/prediction-agreement-partial-v1 --output results/experimental_structures/prediction-agreement-summary-partial-v1
```

Resource planning reserved one CPU thread, 8 GB memory and 1 GB output with a
0.05–2 hour estimate. No new services were used. Experimental quality, local
reliability, complex/conformational context, refinement starting models and
training/template overlap remain unresolved. These are descriptive agreement
measurements, not an unbiased accuracy benchmark or evolutionary distances.

The idle CA mapper now honors explicit excluded entries pinned in the coordinate
retrieval configuration, permitting the full experimental-only expansion while
retaining its four integrative deferrals. Unknown deferrals and mismatched entry
universes are rejected. The two existing CA classification tests still pass.


## Full coordinate retrieval completed and taxon coverage assessed

All 1,032 experimental coordinate archives are complete: 215 verified reused
archives and 817 new downloads, totaling 1,952,380,079 compressed bytes and
8,484,682,674 bytes after decompression. Independent readback checked every
archive checksum, entry receipt, configured entry universe and reuse count.
The producer also streamed full gzip validation. Full CA mapping is now running
in `results/experimental_structures/ca-mapping-full-v1`, retaining all deposited
models, chains, missing positions and ambiguous observations. Planning reserves
one CPU worker, 8 GB memory and 10 GB output with a broad 0.5–8 hour forecast.

```bash
python scripts/map_experimental_ca_residues.py --screen results/experimental_structures/sequence-screen-full-v1 --coordinates data/experimental_structures/coordinates-exact-full-v1 --output results/experimental_structures/ca-mapping-full-v1
```

Reference coverage is concentrated. The 80 exact-sequence candidate proteins
project to 98 marker/taxon links among **nine project taxa, all Ascomycota**.
The remaining 517 analysis taxa have no exact-sequence link in this snapshot.
Saccharomyces cerevisiae accounts for 78 proteins and 1,021 candidate entries;
Aspergillus fumigatus and Schizosaccharomyces pombe each account for one protein.
The six other covered taxa are Saccharomyces entries, including the two curated
hybrids, linked through identical protein sequences. These nine taxa are not
nine independent experimental sources or nine verified unique species.

The coverage table preserves all 526 analysis taxa and their zero-coverage rows.
Taxon, protein and marker counts were independently read back against the 98
exact-sequence links. Deposited organism identity is not inferred from target
sequence identity. No absence of experimental homologs is implied by absence
from this accession-based exact-match inventory. Experimental benchmarks cannot
currently support broad fungal-lineage generalization without further coverage.

```bash
python scripts/summarize_experimental_taxon_coverage.py --references results/experimental_structures/reference-metadata-v2 --snapshot results/structural_markers/gdm-expanded-v1 --output results/experimental_structures/taxon-coverage-v2
```

An unpublished initial coverage calculation used the draft manifest and mutated
a zero-coverage dictionary while rendering, inflating its covered-taxon total.
It is explicitly marked `INVALIDATED.json` in its local output. The published v2
uses the frozen 526-entry analysis manifest and checks coverage against nonzero
rows and distinct links; no v1 numbers are used downstream.


## Starting-model dependencies and prediction chronology

Reviewed all 1,032 experimental entries against their frozen method-specific
metadata. Two entries explicitly list AlphaFold starting models (8RAM, 8RAP),
five list other computational models (4L9P, 4MBG: PHYRE; 9E2W, 9E2X, 9E2Y:
unspecified Other), 548 report only experimental starting models, and 477 lack
starting-model annotations. The raw lists are preserved alongside review flags.
These are entry-level annotations and cannot yet be attributed to the exact
matched target chains. Missing annotations are not evidence that no prediction
was used, and experimental starting models can have their own dependencies.
The authoritative source fields are the frozen RCSB entry responses and
`pdbx_initial_refinement_model` records, including the entry response for
[8RAM](https://data.rcsb.org/rest/v1/core/entry/8RAM).

The experimental-only exact sequence inventory contains 2,843 entity/model rows
(the full inventory also includes seven rows in four integrative entries).
Of these, 2,315 entries were initially released before the linked prediction's
recorded creation date and 528 after it. Dates are stored for each link.
Prediction creation dates are not training or template cutoffs. Later releases
do not establish absence of earlier homologs, sequences, alternate structures,
or related templates. Every reference retains unresolved training independence.

All 1,032 source starting-model lists and the ordering of all 2,843 date pairs
passed readback, along with exact identities of the seven flagged entries.
This review adds concrete dependency evidence; it does not complete chain-level
methods review or establish an independent accuracy test. The current partial
geometry results remain descriptive and should consume these flags during
subsequent benchmark selection and sensitivity analyses.

```bash
python scripts/review_experimental_model_provenance.py --references results/experimental_structures/reference-metadata-v2 --predictions results/structural_markers/gdm-expanded-v1 --output results/experimental_structures/model-provenance-review-v1
```


## ESMFold controls for the experimental reference set

Prepared all 80 experimental-reference sequences for comparison against the
existing local ESMFold configuration. Forty-five canonical full proteins are
at most 512 residues and absent from all four prior immutable queues; the other
35 exceed the current length limit and are explicitly deferred without
truncation. No candidate was chosen by prediction agreement, resolution or
phenotype. All emitted sequence hashes, length/alphabet bounds and disjointness
from original, follow-on, predictor-control and ecology queues passed readback.

Observed timings from 6,666 completed original-run prediction receipts yield
0.078 GPU hours at within-length-bin medians, 0.056–0.106 hours using observed
10th–90th percentile timings, and a planning allowance of 0.158 hours (about
9.5 minutes) after 50% overhead on the upper projection. These are planning
scenarios, not statistical uncertainty intervals. Reserve 24 GB VRAM and 20 GB
disk headroom on the existing host; no new paid services are used.

The experimental-control controller is now waiting for original GPU1 process
691140 with its pinned process start time to finish. It requires that run's
matching configuration and clean completion receipt, zero remaining eligible
sequences, no interruption or OOM, and an idle GPU before launching the 45-control
queue. Inputs, checkpoint, producer and resource receipt are pinned. A launch
receipt prevents duplicate scheduling; failures require explicit review. This
uses a separate copy of the currently live ecology controller with purpose
labels changed; its control logic is identical and the two existing completion/
process-identity tests pass. The live ecology controller remains unchanged.

```bash
python scripts/prepare_experimental_prediction_controls.py --output data/prediction_inputs/experimental-controls-v1
python scripts/advance_experimental_control_predictions.py --config metadata/experimental_control_prediction_controller_config.json
```

The controller is already running; do not start a duplicate. Output will be
`results/predictions/esmfold-experimental-controls-v1`. Independent artifact
validation and matched AlphaFold/ESMFold/experimental comparisons remain pending.
Alternate prediction does not establish experimental or training independence,
and this reference set remains taxonomically narrow.


## Diagnostic figure: per-protein agreement and reference attrition

![Initial experimental agreement](figures/experimental_agreement.svg)

The figure shows every protein with accepted baseline coverage, sorted by its
median CA RMSD, with its baseline number of chain/model comparisons. Missing
higher-threshold points indicate failed coverage, not zero deviation. The
right panel compares the median across all observations with the median across
within-protein medians. Coverage declines from 33 proteins at baseline to 32 at
pLDDT 70 and 18 at pLDDT 90. One of the 34 reference proteins has no eligible
chain/model comparison even at baseline. All threshold row counts and plotted
protein medians were recomputed directly from comparison tables.

These are dependent, partial-snapshot descriptive measurements. Different
thresholds change residues and sometimes the available chain/model and protein
cohorts; fixed-cohort summaries are separately provided above. No visual trend
establishes prediction accuracy or inherited structural change. Experimental
quality, complex context and training independence still require review.

```bash
python scripts/plot_experimental_agreement.py --comparisons results/experimental_structures/prediction-agreement-partial-v1 --summary results/experimental_structures/prediction-agreement-summary-partial-v1 --output results/experimental_structures/prediction-agreement-figure-partial-v1
```

SVG, PNG and PDF exports and a hash receipt are retained in the immutable figure
output. The SVG is published in the repository; the rendered PNG was visually
inspected for readable labels, unclipped axes/legend and correct missing-point
representation.

## Full experimental CA mapping completed and audited

All 1,032 experimental entries completed mapping, retaining 1,306,701 complete
sequence-position records across deposited models and selected chains:

| CA correspondence class | Position records |
| --- | ---: |
| Unambiguous, full occupancy | 1,008,964 |
| No CA observation | 297,499 |
| Multiple CA records | 107 |
| Nonstandard or mismatching monomer | 57 |
| Alternate location or partial occupancy | 73 |
| Invalid coordinates or occupancy | 1 |

Every exported sequence grid, sequence hash, embedded atom identity and count
partition passed full readback. The 1,009,309 embedded atom records include
multiple/alternate observations; their count is not the accepted-residue count.
Independent raw mmCIF field readback covered 2,842 position rows in five
identity-hash-selected entries (7KEE, 3NZX, 7E8S, 5N5Y, 7E93), not every raw
source file. The one invalid-occupancy record was additionally checked directly
in 5LQW: entity 14, model 1, label chain N, residue 640 has occupancy 0.00.
Its coordinates are finite, but it is excluded as an observed CA by policy;
this is not a claim of corrupted deposition.

```bash
python scripts/audit_experimental_ca_mapping.py --mapping results/experimental_structures/ca-mapping-full-v1 --screen results/experimental_structures/sequence-screen-full-v1 --coordinates data/experimental_structures/coordinates-exact-full-v1 --output results/experimental_structures/ca-audit-full-v1
```

The full prediction–experiment comparison is now running against this audited
mapping, using the same fixed protocol as the original partial comparison:
retain all chains/models; predicted focal pLDDT 0/70/90; at least 50 matched CA
and half the full canonical sequence; independently recompute every accepted
geometry with the second rotation/distance implementation. Resource planning
reserves one CPU thread, 8 GB memory and 5 GB output with a 0.1–4 hour forecast
on the existing host. Experimental quality, target-chain context and training
independence remain unresolved; source coverage is still concentrated in the
same nine Ascomycota project taxa.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/compare_experimental_predictions.py --mapping results/experimental_structures/ca-mapping-full-v1 --screen results/experimental_structures/sequence-screen-full-v1 --predictions results/structural_markers/gdm-expanded-v1 --references results/experimental_structures/reference-metadata-v2 --output results/experimental_structures/prediction-agreement-full-v1
```
