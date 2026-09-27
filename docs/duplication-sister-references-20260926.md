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
