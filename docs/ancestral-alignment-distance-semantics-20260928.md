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

## Full frozen-input geometry check

`scripts/ancestral_extant_alignment_geometry.py` validates sampled extant
sequences against observed inputs, retaining an input X as one residue
position even when its sampled amino acid changes. It exports A/gap strings
solely as a positional encoding for alignment comparisons. Internal-only
columns are removed. These strings must never become sequence-inference or
structure-prediction inputs.

An independent signature records, for each extant residue and target tip,
the aligned target residue position or zero for a gap. Twice the number of
changed entries equals the weighted symmetric feature difference above.
Rows follow sorted tip identifiers and ungapped input positions. Comparisons
require the same observed inputs; equal array dimensions alone do not
establish compatible identities. Four tests cover all 733 original fixture
comparisons, input-X resampling, internal-only columns and invalid inputs.

The full check uses all 135 frozen effective input groups and samples at
iterations 0, 10 and 20, with up to 622 extant tips. Resource inventory found
676,652,745 residue-by-tip cells summed across groups; three simultaneous
largest signatures would occupy 2.13 GB. The implementation retains only two
such signatures and is limited to one CPU and 8 GiB memory, without swap.
Runtime planning is 5–60 minutes and output allowance 10 MiB. Launch evidence
is in `metadata/ancestral_extant_geometry_fixture_launch_20260928.json`.

The first attempt failed an exact integer comparison: native output printed
128,405,000 for an independently calculated 128,405,186. Its terminal log and
script are preserved under `results/ancestral/extant-geometry-fixtures-20260928-v1`.
The subsequent check explicitly compares six-significant-digit representations,
retaining exact independently calculated integers, rounded expectations and
raw native stdout. Agreement therefore establishes agreement at the printed
precision, not native bit-exact integer equality. This sensitivity must remain
visible in reporting. Output is under
`results/ancestral/extant-geometry-fixtures-20260928-v2`; require a successful
terminal service and a complete checked receipt before claiming all inputs pass.
