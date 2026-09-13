# Connecting structures to phylogenetic sites

Run `python scripts/map_marker_structures.py --output results/structural_markers/snapshot-v1` to make an immutable snapshot from currently verified AFDB models and the complete 125-marker panel. Choose a new output directory as acquisition advances.

Each linked model must have the exact complete amino-acid sequence of a marker protein. The script rechecks its CIF checksum, polymer sequence and alpha-carbon coverage. Identical sequences can reuse a model across taxa, with species and original protein identifiers retained; shared coordinates are not independent experimental observations. Where multiple candidates have identical sequences, the snapshot chooses the highest mean full-protein pLDDT, then version and model ID as deterministic tie-breaks. Candidate count and selected-model provenance are recorded. This rule is a baseline, not proof of superior biological accuracy.

The full Stockholm alignment connects each original protein residue to its BUSCO profile position. Retained profile positions then map into the concatenated species matrix. Inserted residues advance protein coordinates; gaps do not. The compressed residue table records marker, taxon, protein, model/version, matrix column, protein residue number and CA pLDDT. This makes later site-level sequence/structure comparisons traceable without treating raw structural distances as additive evolutionary branch lengths.

The first executed snapshot screened all 59,840 marker proteins and linked 151 proteins in two taxa, with 57,413 matrix-to-structure residue links. No marker had four linked taxa, so this snapshot does not support branch-specific structural inference. Confidence summaries describe the retained marker residues. PAE, domain boundaries, prediction-source effects, orthology and phylogenetic uncertainty still require assessment.

## Retrieval ordering

Current acquisition first visits marker candidates in round-robin taxon order, then processes every other exact-sequence candidate. No atlas candidates are discarded. The queue snapshot records its size and taxon coverage. Two workers and four outstanding futures keep concurrent requests bounded. New stop handling finishes the few active requests and records their receipts before exiting. Completed CIF checksums are verified on resume.

The earlier queue processed complete proteomes consecutively, delaying comparative coverage. It was intentionally terminated after preserving 6,924 valid model receipt records, then replaced only after its process was confirmed terminal. The new queue started with 12,007 priority marker accessions across 113 taxa, followed by the remaining candidates in a 1,043,924-accession snapshot. Matching continues across the full cohort, so later snapshots will add newly available candidates. The retrieval rate and full-data objective are unchanged.
