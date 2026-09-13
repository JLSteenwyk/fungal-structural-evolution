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
