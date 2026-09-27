# Candidate sequence context — September 27, 2026

The complete 947-candidate annotation set now distinguishes three different
notions of sequence identity: letters paired in each reference-based common
core, the complete unmasked domain interval strings, and the full predicted
protein sequence hashes. All 1,777 interval triads were retained.

| Identity statement | Candidate combinations |
| --- | ---: |
| All aligned common cores have identical duplicate-pair letters | 54 |
| All complete interval pairs have identical sequence strings | 48 |
| Both statements hold | 40 |
| Identical full-protein sequence hashes | 0 |
| Identical complete intervals, but not every aligned core has identical letters | 8 |

These are event/domain combinations, not independent gene pairs or origins.
Among the 54 identical-core cases, 34 are concordant, 15 discordant (14 fungal,
one outgroup) and five sequence-unresolved. Of the 14 fungal discordant cases,
ten also have identical complete intervals and four do not. Conversely, three
additional fungal discordant cases have identical complete interval strings but
not identical paired letters in every common-core alternative.

Complete-interval equality compares exact strings without realignment. Interval
inequality can include different boundaries or lengths as well as substitutions.
Full-chain sequence hashes differ in every candidate, including those with
identical domains; surrounding sequence, domain placement and prediction
variation therefore remain possible explanations. These data do not demonstrate
sequence-independent structural evolution or establish a causal context effect.

The eight identical-interval/nonidentical-core cases expose dependence on
residue correspondence in the reference-based mapping. Their domain strings
match, but some compared residue pairs carry different letters. They require a
control using sequence-position correspondence before interpreting structural
asymmetry. They remain in the ledger and are not silently dropped or relabeled
as erroneous predictions. Even a sequence-position control would not establish
which predicted conformation is correct.

## Reproduction and verification

Run `scripts/summarize_candidate_sequence_context.py --plan
metadata/duplication_candidate_sequence_context_plan_20260927.json`, followed by
`scripts/check_candidate_sequence_context.py --plan
metadata/duplication_candidate_sequence_context_plan_20260927.json --output
<new-readback-path.json>`.

Output directory
`results/structural_comparisons/duplication-candidate-sequence-context-20260927-v1`
contains `candidate_sequence_context.tsv`, `interval_pair_context.tsv` and a
hash-bound receipt. Every source candidate field is preserved, with model IDs,
versions, complete protein lengths and sequence hashes added. All source interval
identities, bounds and sequence hashes remain available in the triad table.

The independent checker read all **3,426** selected full-mask PDB files, checked
hashes, decoded residue letters and checked position lists against the input
manifest. It reconstructed every interval equality, candidate classification
and summary stratum. Full-chain equality uses the already validated source
sequence hashes; full chains were not reparsed in this stage. See
[completion record](../metadata/duplication_candidate_sequence_context_completed_20260927.json).
Resources were one CPU, a 4 GiB memory and 1 GiB output allowance, using existing
models only. GPU prediction remains paused. Position-correspondence controls,
model-source checks and phylogenetic/statistical interpretation remain open.
