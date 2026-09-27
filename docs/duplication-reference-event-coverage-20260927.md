# Reference comparison coverage across duplication events

Every one of the 73,888 source event/reference/duplicate-side links now has
explicit availability under both full and pLDDT70 masks: 147,776 link/mask rows.
The two duplicate sides are joined into 73,888 event/reference/mask rows, with
all tied nearest references and both guide alternatives retained.

| Guide | Mask | Event/reference cases | Both sides, both orders numerically usable |
| --- | --- | ---: | ---: |
| Profile | Full | 18,478 | 18,156 |
| Profile | pLDDT70 | 18,478 | 16,673 |
| MAFFT | Full | 18,466 | 18,143 |
| MAFFT | pLDDT70 | 18,466 | 16,653 |

These are availability counts, not independent events, tested effects or
validated asymmetry cases. Multiple tied references, masks and guides reuse
underlying genes and models. The source inventory comprises provisionally
eligible references only; original events lacking a provisional reference remain
explicit in the upstream complete sister-reference inventory.

The original source links comprise 72,598 additional-pair comparisons, 641
comparisons reusing a primary duplicate pair, and 649 identical-model links.
The additional pairs use the completed audited numerical summaries. Reused
primary pairs remain `awaiting_primary_pair_audit`; the current projection does
not ingest incomplete checkpoints. Identical-model links remain
`identical_model_not_independent_comparison` without invented structural scores.
All source identities, the deterministic representative flag and alternative
tied references remain in the expanded table. No native direction is selected.

## Reproduction and verification

Run `scripts/project_reference_event_coverage.py`, then
`scripts/readback_reference_event_coverage.py`. The output directory is
`results/structural_comparisons/duplication-reference-event-coverage-20260927-v1`.
Preserve existing outputs and use new paths for reruns.

The independent dataframe checker reconstructed every original field, both-mask
expansion, model-pair availability join, duplicate-side pivot and all 24 status
combination summary rows. Full proof and summary are versioned as
`metadata/duplication_reference_event_coverage_{receipt,readback}_20260927.json`
and `metadata/duplication_reference_event_coverage_summary_20260927.tsv`.

Numerical availability still requires coverage/confidence qualification and
shared-residue comparisons before interpreting asymmetry. Extant sister
references are not reconstructed ancestors or established orthologous outgroups.
