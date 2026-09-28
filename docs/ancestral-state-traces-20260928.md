# Residue-anchored ancestral-state traces

The ancestral uncertainty workflow now preserves categorical state traces,
in addition to scalar logs, candidate lengths and extant alignment distances.
For each saved alignment and candidate ancestor, each extant input residue
anchors the column from which the ancestor's amino-acid or gap state is read.
Coordinates are tip identifier plus one-based position in the ungapped input
sequence; domain inputs are not automatically whole-protein coordinates.

This preserves changes in amino-acid state and alignment assignment without
assuming fixed alignment columns across samples. Multiple extant residues
can anchor the same ancestral residue in one sample. They are correlated
views of the same reconstruction, not independent observations or sites.
The alphabet explicitly includes X and gap. A gap at an anchor is not by
itself evidence of loss of the whole protein or domain.

For each chain, `states.npz` contains the complete saved-iteration × source-node
× extant-residue state array, iteration numbers, unanchored ancestral residue
counts and state-count arrays after discarding iterations through 250 and
500. The initial horizon provides 75 and 50 retained samples respectively.
`coordinates.json` records alphabet order, stable source-node identifiers,
tip order, sequence lengths and the hashed biological input alignment.
Full temporal traces remain available for autocorrelation and indicator-based
diagnostics. Counts alone cannot supply effective sample sizes.

Ancestral residues in columns without extant anchors have no established
cross-sample identity. Their counts are retained separately; the original
sampled alignments preserve their complete sequences and ordering. No anchor
is invented for these columns.

Projection is checked independently against direct indexing of the ancestral
aligned sequence at every extant tip's nongap columns. Candidate strings are
also reconstructed, unanchored residue counts independently checked, and
all arrays read back after serialization. State-count totals must equal the
number of retained samples at every node/anchor for both cutoffs. Three
synthetic tests check shared anchors and gaps, an unanchored insertion and
changed runtime node label, original-X resampling, and rejection of a changed
known extant residue.

The first application uses the frozen seven completed-chain snapshot. It
retains all 707 saved alignments; it does not select only stable-looking
regions. The resource plan limits the run to one CPU, 4 GiB memory and no
swap, with 1–30 minutes runtime and 1 GiB output allowances. Production
chains and pinned scripts are unchanged.

```bash
python scripts/prepare_ancestral_state_traces.py \
  --snapshot metadata/baliphy_first_completed_chains_20260928.json \
  --inputs results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json \
  --output results/ancestral/first-anchored-state-traces-NEW
```

Results are under `results/ancestral/first-anchored-state-traces-20260928-v1`.
These arrays are diagnostic inputs, not qualified posterior probabilities.
Four-chain comparisons, rare/constant-state handling, alignment mixing,
adequacy and topology/root sensitivity remain required before ancestral
structure ensembles can be qualified. One reference-anchor projection cannot
capture every mode of the joint alignment/ancestral-sequence posterior.

The first application completed in 46.6 seconds. All **47,654,628 recorded
state observations** passed the independent indexing comparison; both
burn-in count arrays and full temporal arrays passed serialization checks.
The service was inactive/success/exit 0, and artifact and input/script hashes
were verified after termination. Completion evidence is
`metadata/baliphy_first_anchored_states_completed_20260928.json`.
These counts include correlated repeated views of ancestral residues, four
candidate ancestors and all saved iterations; they are not a count of
independent biological sites or posterior samples.
