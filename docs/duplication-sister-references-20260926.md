# Provisional sister-clade references for duplicate-pair comparisons

This full-cohort inventory supports the planned structural-asymmetry work.
Both guides' complete reviewed terminal-candidate sets are included: 109,245
profile and 109,228 MAFFT event rows. It does not alter the running frozen
duplicate-pair alignment queue or infer new structures.

For each exact duplicate tip pair, the procedure inspects the immediate
parent and all sister-clade descendants. Every candidate retains sister-gene
and taxon counts, focal-taxon presence, model availability, parent branching
and whether that parent appears in the native duplication table. A
provisional reference is available only if the parent is bifurcating, is not
reported as a duplication, the sister clade has no focal-taxon genes, and
at least one nonfocal sister protein has a frozen structure model.

Missing parents, reported parent duplications, multifurcations, focal-taxon
sister genes and absent models remain explicit dispositions. The absence of
a reported parent duplication does not establish speciation. These are
extant reference candidates, not validated orthologous outgroups or
reconstructed ancestral states. Hidden paralogy, differential loss, uncertain
rooting and topology still need review.

The nearest modeled nonfocal sister tip is identified using the saved
sequence-tree path length from the duplicate node. All ties within absolute
1e-12 are retained; the lexical first tied gene is a deterministic provisional
representative. Selection uses sequence-tree distance and model availability,
not observed structural divergence. Path lengths to both duplicate tips and
the duplicate-pair path are retained as relative sequence divergence, not
calendar time. References are recorded even when another eligibility flag
prevents provisional use, so these cases remain inspectable.

Implementation: `scripts/inventory_duplication_sister_references.py`.
Tests: `python scripts/check_duplication_sister_references.py`. Known-tree
checks passed expected reference eligibility and four exact path distances;
tied tips, zero branch lengths, missing parent, multifurcation, reported
parent duplication and focal-taxon sister cases were checked. Negative
branches and duplicate labels were rejected.

Plan: `metadata/duplication_sister_reference_plan_20260926.json`.
Service: `fungal-duplication-sister-references-20260926.service`.
Launch identity: `metadata/duplication_sister_reference_launch_20260926.json`.
Output: `results/orthology/duplication-sister-reference-inventory-20260926-v1/`.
The full scan uses one CPU, 8 GiB RAM and no swap, with a 1-GiB output
allowance and uncalibrated 0.1–8-hour planning range. Source event tables,
resolved trees, model bridge and completed candidate review are pinned and
rechecked. All candidate rows must receive dispositions before completion.

The inventory is running. Independent readback, cross-guide reference
agreement, reference model quality and domain checks, structural comparisons
to references, alternative-reference sensitivity and asymmetry tests remain
pending. This step does not establish any duplication-associated effect.


## Completed inventory and guide comparison

Both full scans completed and the service exited successfully. Every candidate
received a disposition:

| Immediate-reference disposition | Profile | MAFFT |
| --- | ---: | ---: |
| Provisional reference available | 17,461 | 17,448 |
| Parent reported as a duplication | 54,891 | 54,884 |
| No modeled nonfocal sister protein | 36,817 | 36,820 |
| Parent not bifurcating | 76 | 76 |
| Total candidate rows | 109,245 | 109,228 |

The categories follow the documented priority order; for example, an event
with a reported parent duplication may also lack a modeled reference. They
are not mutually exclusive biological explanations. The table reports the
single output disposition for each candidate.

`scripts/compare_duplication_sister_references.py` independently checked all
original candidate columns and the complete protein-pair universe against
the prior reviewed exports. Among 108,918 shared duplicate pairs, 17,392 are
provisionally eligible under both guides. For 17,365, the complete nearest
reference sets agree and the chosen reference gene/model also agrees. The
remaining 27 have disjoint nearest-reference sets; none have partially
overlapping sets. These cases remain explicit for later reference sensitivity.

Among shared pairs, 28 are provisionally eligible only in profile (13 lack a
modeled reference in MAFFT; 15 have a reported parent duplication there).
Nineteen are provisionally eligible only in MAFFT (one lacks a modeled
reference in profile; 18 have a reported parent duplication there). Separately,
41 profile-only and 37 MAFFT-only duplicate candidates have provisional
references. The complete eligibility matrix is
[`metadata/duplication_reference_guide_eligibility_matrix_20260926.tsv`](../metadata/duplication_reference_guide_eligibility_matrix_20260926.tsv).

Source-plan pins, receipts, all output artifact hashes and count identities
were rechecked and archived in
`metadata/duplication_sister_reference_completed_20260926.json`. This comparison
does not independently reconstruct sister-clade selection or path distances.
High cross-guide agreement is not proof of orthologous reference assignment
or a structural-asymmetry effect. Most candidate pairs still need deeper
gene-tree/context review or additional reference coverage; they remain in
the broader duplication project.

Reference-set fixtures passed shared/guide-only, unavailable, identical,
overlapping and disjoint cases; an empty provisional set was rejected.
Reproduce with `python scripts/check_duplication_reference_comparison.py`.
Reproduce the full comparison with:

```bash
python scripts/compare_duplication_sister_references.py \
  --inventory results/orthology/duplication-sister-reference-inventory-20260926-v1 \
  --plan metadata/duplication_sister_reference_plan_20260926.json \
  --review results/orthology/duplication-candidate-tree-review-20260926-v1 \
  --output <fresh-output-directory>
```

Current comparison output is
`results/orthology/duplication-reference-guide-comparison-20260926-v1/`.
