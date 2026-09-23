# Whole-proteome structural domain interval registry

The registry links every one of the frozen catalog's 1,290,278 AlphaFold models
to the independently audited full Pfam candidate-architecture database by exact
sequence hash. It retains the 1,319,513 original species/protein/model links,
checking each against the annotation protein table as well as the model table.
No prediction, structural clustering or evolutionary inference is performed
by this stage.

`structure_domains.sqlite` contains four related tables:

- `models`: structure identifier/version, sequence hash, length, path, raw-hit
  count and policy-agreement flag, including models without hits.
- `segments`: every distinct retained annotation hit across policies, including
  Pfam accession/version, type, clan, alignment bounds, envelope bounds and HMM
  coverage. Repeated hits remain separate through their source hit identifiers.
- `policy_hits`: retained hit membership for all four competition policies,
  conservative-architecture flags and domain-interval candidate flags.
- `protein_links`: every original taxon/protein association to its model.

Coordinates are **one-based inclusive source-protein positions**, with alignment
and envelope boundaries retained separately. Each interval must satisfy
`1 <= envelope_start <= alignment_start <= alignment_end <= envelope_end <= protein_length`.
Inconsistent same-hit fields across policies, repeated hits within a policy,
unknown model sequences, mismatched protein sequences and invalid bounds stop
the build. Database keys preserve distinct repeated occurrences.

Conservative architecture requires policy agreement and absence of rank ties,
unresolved overlaps, candidate nesting, alignment overlaps and retained partial
HMM hits below 0.70 coverage. A domain-interval candidate additionally requires
Pfam type `Domain`, HMM coverage >=0.70 and alignment span >=30 residues. These
are screening choices, not validated structural boundaries. All other types,
intervals and policy alternatives remain available for sensitivity analyses.
No-hit models remain explicitly unannotated; they are not domain-loss calls.

This registry does not extract coordinates or apply residue pLDDT/PAE criteria.
Those checks, independent registry readback, domain-level structural comparisons,
clustering and phylogenetic interpretation remain required. Exact sequence
identity supplies the intended residue numbering; it does not independently
validate every CIF coordinate here.

## Execution and resources

Script: `scripts/build_structure_domain_registry.py`.
Plan: `metadata/whole_proteome_structure_domain_registry_plan.json`.
Output: `results/domains/whole-proteome-structure-domain-registry-20260923-v1`.
Unit: `fungal-structure-domain-registry-20260923.service`.

The plan pins source data and code. Its catalog and architecture checksums were
also matched to the previously completed catalog and family-domain bridge
receipts before launch. Limits are one CPU, 16 GiB memory and no swap. Planning
allows 10 GiB output, at least 100 GiB free disk and an uncalibrated 0.5–12 hours;
no GPU or paid resources are used. The process was confirmed live in
`metadata/whole_proteome_structure_domain_registry_launch.json`.

Fixtures checked union-hit deduplication with preservation of four policy
memberships, conservative exclusion on policy disagreement and rejection of an
out-of-bounds envelope. A completed producer receipt will still be labeled
pending independent readback. At this checkpoint the build is running and no
registry result is claimed.
