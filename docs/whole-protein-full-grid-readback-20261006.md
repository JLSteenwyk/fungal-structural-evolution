# Full whole-protein numerical readback

The independent full-grid readback completed on 2026-10-06. It checked every
one of the 1,673 serialized numerical-followup results drawn from the complete
375,350-disposition whole-protein grid.

The final receipt is
[`whole-protein-full-flag-readback-20261006-v2/receipt.json`](../results/model_validation/whole-protein-full-flag-readback-20261006-v2/receipt.json).
It records 1,673 unique `(fit_input_id, tree)` identities and 1,673
`independently_replayed_numerical_followup` outcomes. The readback manifest
has SHA-256 `5c30c6ac2e9d47b8824d01aa7ebb90688e28e69daf5b60aa67ed1aa963095e38`.

The superseded reader stopped because a fixed central finite-difference step
of `1e-5` was truncation-sensitive for one interior optimum. The replacement
reader used the predeclared `adaptive_v2` policy: central differences at
`1e-6`, while retaining every source optimisation-review disposition. The
maximum independently reconstructed objective discrepancy was
`7.639755494892597e-11`.

This validates the serialized follow-up arithmetic, candidate likelihoods,
score checks, coefficient/covariance transforms and recovery lineage. It does
not establish a global likelihood optimum, calibration, structural homology,
or a biological sequence--structure conclusion. Those gates remain separate.

## Complete fit/audit/follow-up closure

The receipt-based closure subsequently bound the complete 75,070-input,
375,350-disposition fit and audit census plus every 1,673 follow-up. Its
[completion receipt](../metadata/full_whole_protein_optimization_completed_20261006_v2.json)
records 527,371 bound source hashes, 41,835 independently replayed candidate
likelihoods, and a maximum objective reconstruction error of
`7.639755494892597e-11`. The hash archive is stored outside Git at
`results/model_validation/full-whole-protein-optimization-closure-20261006-v2/`.
Three original producer-service journals were verified; the collected v2
reader is established by its immutable receipt and hash-checked manifest.
