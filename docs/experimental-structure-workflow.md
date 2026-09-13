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
