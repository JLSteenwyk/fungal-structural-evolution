# Within-family candidate annotation variation

This full-data inventory describes candidate annotation differences within
both expanded family partitions. It prepares family-level inputs for domain
evolution analyses; it does not infer evolutionary events or their direction.
The execution waits for the complete independent family/domain bridge audit.

Each family receives four rows, one for each alignment/envelope and
E-value/bit-score competition policy. Rows retain protein, taxon and unique
sequence counts; taxa with multiple family members; no-hit counts; policy
disagreement; rank ties; unresolved overlaps; candidate nesting; alignment
overlap; and partial HMM matches. Multiple family members in a taxon are not
automatically called duplications: reconciliation and annotation review remain
necessary.

For proteins with retained annotations, the inventory separately counts:

- Ordered model/type signatures, preserving repeated annotations.
- Unordered multisets, preserving annotation multiplicity.
- Model/type sets, ignoring order and multiplicity.
- Multisets with more than one observed ordering.
- Model/type sets with more than one observed multiplicity pattern.
- Taxa containing more than one signature within a family.

These categories distinguish descriptive composition, multiplicity and order
variation. They are not estimates of fusion, loss, duplication or rearrangement
events. All Pfam annotation types remain explicit; a Family or Repeat model
is not relabeled as a Domain. No-hit proteins are counted as missing annotation
evidence, never as an empty architecture supporting a loss.

Each row also contains a conservative sensitivity summary excluding proteins
with policy disagreement, rank ties, overlaps, candidate nesting, or any
retained HMM match with coverage below 0.70. This is a declared screening
choice, not proof of correctness. It may exclude real nested or partial
architectures and must be compared with the complete observed summary.

Source family indices differ between guides. Compare exact memberships using
the bridge `families.membership_sha256` and source crosswalk; do not equate
matching family names. Multi-copy taxa, incomplete annotations and uncertain
gene homology need further checks before biological interpretation.

## Execution and validation

`scripts/summarize_family_architectures.py` reads all 5,815,847 protein links
under each guide, without reducing taxa or families. It streams one family at
a time and writes compressed TSVs under
`results/domains/family-architecture-variation-v1/`. The final receipt records
source hashes, full family/protein counts and output checksums.

The resource plan is `metadata/family_architecture_variation_plan.json`:
one CPU equivalent, 16 GiB RAM, no swap, 10 GiB output allowance, and a
100 GiB free-disk gate. The 1–24 hour allowance is uncalibrated planning,
not a measured ETA. The service is
`fungal-family-architecture-variation-20260922.service`; it waits for the
exact bridge-auditor process and requires a passing receipt before analysis.
No GPU prediction or paid resource is used.

`python scripts/check_family_architecture_variation.py` verifies separation of
order and multiplicity variation, taxon-specific copy counts, model-type
preservation, no-hit exclusion, and conservative screening. Full output
readback and guide sensitivity comparisons remain required after execution.
Subsequent event analyses must use supported genealogies and reconciliation,
check sequence/annotation quality, account for missingness and uncertainty,
and distinguish within-domain structural change from architecture turnover.

## Updated dependency after bridge audit performance recovery

The initial waiting job was stopped without running the family analysis while
the bridge audit was replaced to use its existing family index. The unchanged
architecture summarizer already explicitly uses that index. It is queued again
under `fungal-family-architecture-variation-recovery-20260922.service`, with
`metadata/family_architecture_variation_recovery_plan.json` and fresh output
`results/domains/family-architecture-variation-v2/`. This plan binds the new
audit PID/start ticks and requires the successful v2 bridge readback. CPU,
memory, disk and validation requirements are unchanged. The prior incomplete
output directory remains preserved. See the [bridge recovery](family-domain-bridge-20260922.md#indexed-readback-recovery).

## Full summary readback

`scripts/readback_family_architecture_variation.py` is queued behind the exact
replacement producer process under
`fungal-family-architecture-readback-20260922.service`. It requires the completed
producer receipt, the passing independent bridge audit and matching source
database hashes. Its pinned plan is
`metadata/family_architecture_variation_readback_plan.json`; fresh output is
`results/domains/family-architecture-variation-readback-v2/`.

For every family and each of the four policies, the audit independently
reconstructs protein/taxon/sequence counts, missing-hit and uncertainty counts,
ordered signatures, multiplicity-preserving multisets, model/type sets and
within-taxon variation. It implements multiplicity as sorted repeated tokens,
separately from the producer's token-count representation, and does not import
the producer's metric functions. It checks the complete observed and
conservatively screened summaries, exact row order, absence of extra rows and
full partition totals. The SQLite input join is shared in design; this audit
does not independently repeat Pfam searches or infer family membership.

Pre-launch checks passed with `python scripts/check_family_architecture_readback.py`:
known order/multiplicity examples, 100 seeded cases varying repeats, annotation
types and QC flags, rejection of altered/missing/repeated rows, and an end-to-end
two-guide database fixture. These tests validate the implementation, not the
unfinished production output. The resource allowance is one CPU, 16 GiB RAM,
no swap, 1 GiB output and 50 GiB free disk, with an uncalibrated 1–24 hour
planning interval after the producer finishes. The successful output status
will be `passed_full_family_architecture_variation_readback`; no such production
completion is claimed at launch.

## Guide sensitivity with exact family memberships

The bridge readback passed, and the v2 architecture producer is running.
`scripts/compare_family_architecture_guides.py` is queued after the independent
summary readback under
`fungal-family-architecture-guide-comparison-20260922.service`. Its plan is
`metadata/family_architecture_guide_comparison_plan.json`; output will be
`results/domains/family-architecture-guide-comparison-v2/`.

The comparison matches families by the hash of their complete sorted native
gene membership, not by their guide-specific family names. For exact shared
memberships, every policy-specific metric must agree, even when family names
differ. Missing, duplicate or mismatched rows fail the comparison.

Per-guide/per-policy summaries report all families, exact shared memberships
and guide-specific memberships separately. They count missing annotations, QC
flags, ordered-signature variation, multiplicity variation and order variation.
The observed and conservative subsets remain separate. Counts also distinguish
families with at least two annotated taxa; this is descriptive coverage, not
proof of replicated evolutionary transitions. No family or taxon is removed
from the full analysis. Guide-specific membership does not itself indicate a
domain evolutionary event. Summing across the all-family and subset rows, or
across policies, would double-count proteins.

The pre-launch fixture passed exact-membership matching across renamed family
IDs, full partition totals, order/multiplicity summaries and rejection of
altered, missing or duplicated rows. Run
`python scripts/check_family_architecture_guide_comparison.py`. The queued
comparison uses one CPU, 16 GiB RAM, no swap, a 1-GiB output allowance and
a 50-GiB free-disk gate; 0.5–8 hours is an uncalibrated planning interval.
Results remain pending until the full producer and audit finish.
