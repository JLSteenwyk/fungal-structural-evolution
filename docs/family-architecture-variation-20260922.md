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
