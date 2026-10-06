# Gene–repeat proximity tables

This full-panel stage joins parsed RepeatMasker calls to the immutable,
source-bound gene-order coordinates on matching sequence identifiers. It emits
all overlapping calls and every tied nearest upstream and downstream repeat;
it does not select one repeat arbitrarily when several share the same boundary.

The per-gene table retains zero, one, or many overlaps, nearest distances and
tie counts, coordinate status, and an explicit no-repeat-call-or-unresolved
sequence-identifier disposition. The latter does not claim a sequence contains
no repeat because a RepeatMasker `.out` table only lists sequences with calls.

These are coordinate covariates for later phylogenetic models. They do not
establish repeat-mediated rearrangement, transposition near a gene, gene
duplication, or protein structural divergence.

`repeat_proximity_full_panel_20261006_v1.sbatch` is the complete 526-taxon
array definition. Its worker checks the frozen taxon manifest, the completed
RepeatMasker receipt and `.out` checksum, and the full-panel gene-order
readback before parsing or linking a taxon. It is deliberately gated on the
independent full-panel RepeatMasker audit; it must not be submitted while that
annotation stage remains incomplete. Each worker refuses to overwrite a prior
output, retains every parsed call and tied nearest relation, and writes hashes
for all downstream artifacts in its completion receipt.

`scripts/audit_repeat_proximity_full_panel_v1.py` independently verifies the
complete panel before these tables become eligible for phylogenetic modeling.
It requires the immutable gene-order full-panel readback, each manifest taxon,
all five expected outputs and their hashes, source identity, coordinate-table
binding, and consistency of parent versus child row counts. Its aggregate
counts are integrity bookkeeping, not estimates of repeat burden or evidence
for repeat-mediated structural evolution.
