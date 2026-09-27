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
