# Residue coordinates for sampled ancestral alignments

Sampled alignment column numbers cannot identify the same position across
chains: insertions and homology changes shift columns. Runtime internal-node
names also differ between processes. `scripts/ancestral_residue_anchors.py`
therefore represents each column by the sorted set of extant input residues
it contains, using `(tip identifier, one-based ungapped input position)`.
Candidate states use source-node identities supplied by the audited tree
mapping. Input positions in cropped/domain alignments are not automatically
whole-protein coordinates.

An anchor set records the sampled homology hypothesis. A change in that set
must remain visible even if its amino-acid letters are identical. A shared
anchor alone does not prove the same homology relationship for the other
residues. Columns lacking extant anchors remain explicit, ordered records
without invented cross-sample identities; their ancestral states are retained
so that the full candidate sequences can be reconstructed.

This is a coordinate-projection component, not a convergence diagnostic.
Production chains are unchanged. Next, comparisons across the four seeds
need to assess ancestral states and homology relationships at the same extant
residue coordinates, while reporting unanchored ancestral material separately.
Scalar and sequence-length diagnostics alone cannot qualify those posteriors.

Synthetic tests cover column shifts, changed homology among identical
letters, runtime-node renaming, unanchored residues, ambiguous input residues,
and malformed alignments. The real-data fixture check covers all 135
effective short-run inputs and all saved iterations (0, 10, 20), checks every
extant coordinate once in order, and reconstructs all four candidate
sequences. Its evidence is recorded in
`metadata/ancestral_residue_anchor_fixtures_20260928.json`.

```bash
python scripts/test_ancestral_residue_anchors.py
python scripts/check_ancestral_residue_anchor_fixtures.py \
  --output metadata/ancestral_residue_anchor_fixtures_20260928.json
```

The fixture command requires a new output path and the checksum-bound local
short-run datasets. These checks validate coordinate extraction; they do not
establish stationarity, adequate effective sample size, or convergence.
