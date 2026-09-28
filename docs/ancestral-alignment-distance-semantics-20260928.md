# Alignment-distance qualification for ancestral sampling

The installed BAli-Phy 4.3 `alignment-distances AxA --distances pairwise`
matched an independently expressed residue-coordinate feature calculation
for all 733 matrix entries in two deterministic fixture groups. This
qualifies the canonical-amino-acid/gap fixture semantics, not posterior
convergence or the full production input alphabet.

For each unordered pair of extant tips, represent aligned residues by their
ungapped sequence positions. Retain residue–residue features with weight 2
and residue–gap features with weight 1. Sum the weights of features present
in only one alignment. Thus `AA/A-` versus `AA/-A` has distance 6: two
changed residue pairings contribute 4 and two changed residue–gap features
contribute 2. Repeated amino-acid letters do not identify homologous residues.
The metric is unnormalized and depends on sequence lengths and tip count.
Do not compare its raw magnitudes across families as equivalent scales.

The fixtures cover randomly generated alignments of identical ungapped
sequences, duplicate samples, a redundant all-gap column and changed FASTA
record order. The independent implementation rejects ambiguous characters
and unequal aligned lengths explicitly. Production input X residues can
have different sampled amino acids: their handling must be qualified using
stable observed residue coordinates before applying the installed tool.
Do not substitute changing ancestral sequences into a metric requiring
identical ungapped sequences across samples.

The source audit binds the installed binary and three downloaded source
files by SHA-256. In the
[version-matched source](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/src/alignment/alignment-util.cc),
`pairs_distance` sums both directions of `asymmetric_pairs_distance`;
`A_match` tests residue coordinate or gap agreement. The `AxA` dispatcher
uses this function. The separately named `NxN` routines are not an
interchangeable definition of this statistic.

Reproduce with a new receipt path:

```bash
python scripts/check_baliphy_alignment_distance_semantics.py \
  --output metadata/baliphy_alignment_distance_semantics_20260928.json
```

The receipt preserves every fixture, expected and observed matrix, source
commit, binary/source/script hashes and runtime. The check took under one
second locally; it does not estimate the cost of full posterior comparisons.

Next requirements are full extant projection and ambiguity validation,
within-chain versus between-chain comparisons for each independently seeded
quartet, and diagnostics of ancestral states conditional on stable residue
anchors. Account for autocorrelation and compare both prespecified burn-in
cutoffs. Alignment distances, scalar diagnostics and length diagnostics are
complementary evidence; agreement in one does not establish convergence in
all. No running chain, pinned script, or production plan was modified.
