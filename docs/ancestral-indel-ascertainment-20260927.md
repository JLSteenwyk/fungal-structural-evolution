# Indel observability and ascertainment

All 156 SIC encodings were inspected by
`scripts/audit_ancestral_indel_ascertainment.py`. Results and per-character
counts are in `results/ancestral/indel-ascertainment-audit-20260927-v1`.
Input receipts, matrices and inspected software sources are checksum-bound.

Of 12,957 characters, 12,117 contain at least one unknown state. Every
nonempty character retains at least one observed 1. Three encodings have no
characters: OG0000294 alignment/FAMSA, alignment/MAFFT and envelope/FAMSA,
all under the terminal-unknown policy. These are explicit empty-data cases,
not evidence for zero historical indel events. No taxa or characters were
removed in this audit.

## Source inspection

In FastML 3.11, `sequenceContainer::startZeroSequenceContainerGL` constructs
one fully observed all-zero pattern when minNumOfOnes=1 and minNumOfZeros=0.
`unObservableData::setLforMissingData` evaluates its probability on the full
tree. `likelihoodComputation.cpp` divides likelihoods by one minus that
probability. This is a common all-taxa ascertainment denominator, not a
character-specific missing-mask denominator.

`gainLoss::checkMinNumOfOnesOrZeros` can lower the minimum and disable the
correction if an input column has no observed 1. Our nonempty matrices do
not trigger that particular condition. This does not resolve the distinction
between conditioning on a gap somewhere in the complete latent pattern
and conditioning on a gap among the observed states.

For a fixed missing mask M, conditioning on an observed 1 would instead use
P(pattern)/(1-P(all observed tips are 0; tips in M unconstrained)). With
multiple rate categories the exclusion probability must be integrated over
categories before conditioning the mixture. Because the selection event is
already implied by any observed pattern containing 1, this normalization
changes parameter fitting but does not further change node posteriors at
fixed parameters. This identity provides a useful implementation check.

However, SIC unknown states arise partly from containing gaps and unknown
flanks; the mask is not demonstrably independent of the latent gap process.
A mask-conditional likelihood is therefore an explicit working model, not
proof that ascertainment is fully solved. Binary SIC characters also overlap
and cannot be treated as independent historical events or independently
sampled into an ancestral sequence without compatibility checks.

## Next inference requirements

Retain the authors' all-taxa correction as a named method sensitivity and
evaluate the observed-mask conditional model separately. Preserve all three
empty-data dispositions. Do not silently replace unknown states by 0, remove
containing gaps, or conflate absence of an exact SIC gap with residue presence.
Compare fitted parameters and node probabilities under these assumptions,
independently check likelihoods and posteriors, and resolve overlapping-gap
compatibility before assembling ancestral sequence ensembles. No indel
model has been fitted by this audit.
