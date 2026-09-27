# Whole-protein duplicate/reference comparison inventory

The full provisional-reference inventory now resolves 36,944 event/reference
links into 17,619 distinct oriented model triples: duplicate A, duplicate B and
an extant reference. All gene identities, guide alternatives, tied references
and deterministic representative flags are retained in the event table.
Computational reuse does not remove species-specific models or turn alternative
links into independent observations.

Of the distinct triples, 15,713 contain three different models, 1,904 contain
two, and two use one model in all roles. Duplicate A/B share a model in 1,585
triples. Reference identity is also explicit (159 A/reference and 166 B/reference
instances, overlapping other identity categories). Model identity is not
independent evidence of structural invariance; no structural scores or residue
mappings have been invented for these cases.

Every A–B, A–reference and B–reference comparison has an exact hashed pair key
and an explicit source: existing primary queue, additional reference queue or
identical model. All nonidentical A–B pairs resolve to the complete primary
queue. The reference paths also retain the cases that reuse primary pairs.
Future evaluation of both confidence masks and all eight combinations of input
order will contain 281,904 dispositions before any eligibility exclusions.

This stage does not fit structures. Its purpose is to establish the complete
input universe for shared-residue comparisons. The next mapping stage must
retain both reference-intersection triples and the subset consistent with the
direct duplicate-pair alignment, carry all exclusions and preserve original
protein positions. Independent complete primary alignment validation remains
a dependency. Identical-model masks require explicit input/geometry handling;
identity does not override insufficient residues or a prior exclusion.

## Reproduction and verification

Run `scripts/prepare_whole_protein_reference_triads.py`, followed by
`scripts/readback_whole_protein_reference_triads.py`. Preserve immutable outputs
in `results/structural_comparisons/whole-protein-reference-triads-20260927-v1`.
The independent checker joined both source-side tables to every event link and
verified all model/version/sequence hashes, oriented triple identities, pair
identities and work partitions. It covered all 73,888 original side links;
no reference was selected or omitted based on structural outcomes.

Evidence: `metadata/whole_protein_reference_triads_{receipt,readback}_20260927.json`.
The original no-reference event dispositions remain in the upstream full sister-
reference inventory. Extant references are not ancestors or proven orthologous
outgroups. No structural asymmetry, selection or ancestral change is inferred.
