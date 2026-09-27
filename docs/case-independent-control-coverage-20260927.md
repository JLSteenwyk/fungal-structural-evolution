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
