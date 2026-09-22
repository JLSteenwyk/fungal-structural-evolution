# Hierarchical orthogroup output verification

The profile-guide native reconciliation has completed, and exact resolved-tree
membership checks passed for 70,307 trees. Hierarchical orthogroup (HOG)
outputs require separate validation before duplication and domain analyses.
The full profile output contains 525 node tables totaling 6,381,538,832 bytes.

`scripts/audit_hog_output_identities.py` is running against every table under
`fungal-profile-hog-identities-20260922.service`. The plan
`metadata/profile_hog_identity_audit_plan.json` pins all tables, native-stage
receipt, source species/protein/family maps, labeled species tree and audit
implementation (532 files). It checks:

- The complete expected set of internal-node tables and source-taxon headers.
- Valid, unique and consecutive node-specific HOG identifiers.
- Every emitted protein's exact taxon and original family identity.
- Containment of every assigned species within the named species-tree clade.
- No repeated gene assignment within one hierarchical level, while permitting
  expected reuse across different levels.
- Stable output files during reading and unchanged source hashes afterward.

The verifier independently reads source membership and does not call the
native HOG inference code. The installed native writer uses all source-species
columns at every level, so containment is tested on nonempty cells rather than
by restricting table headers. Native two-/three-gene HOG output uses a separate
routine that continues past singleton families; the previously observed
small-family orthologue-pair omissions must not be assumed to apply to HOGs.

The output will be `results/orthology/profile-hog-identities-20260922-v1`, with
per-node row/assignment counts and hashes. An assignment count summed across
levels is not a number of distinct proteins. The checks do not establish
complete assignment of all eligible genes, correctness of ancestral group
membership, gene-tree parent-clade labels, hierarchical nesting, rooting,
orthology or duplication events. These remain further requirements, as does
verification of the still-running MAFFT-guide reconciliation.

A hand-constructed three-taxon fixture permits a gene to recur across levels
and rejects wrong-family assignments, genes outside the node's species clade,
within-level duplicate genes, unknown protein identities and a missing table.
Run it with `python scripts/check_hog_output_identities.py`. The fixture passed;
full biological-output verification is still running.

Pre-launch resources are one CPU, 32 GiB RAM, no swap, 0.01 GiB output, and an
uncalibrated 0.1–4 hour planning range on existing local resources. No GPU or
paid service is used. The native output and running reconciliation are read-only.
