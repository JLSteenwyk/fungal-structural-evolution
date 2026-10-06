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
