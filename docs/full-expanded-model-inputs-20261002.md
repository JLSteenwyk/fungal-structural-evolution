# Full expanded sequence–structure model inputs

The complete input stage closed on October 2 for **75,188 logical cases**
and **37,600 physical target–control combinations**, preserving all **4,250,692
original selection memberships**. The [final completion receipt](../metadata/full_expanded_model_inputs_completed_20261002.json)
binds full production, independent readback and provenance closure. This is an
accepted input handoff, not a variance fit or calibrated evolutionary result.
All eight scientific aims remain incomplete.

All 751,880 rows and 51,840 fixed setting counts passed independent readback,
including **10,526,320 numeric-cell checks** against 50-digit Decimal and raw
membership reconstruction in SQLite. The original closure verified **2,369,865
source/artifact bindings** and both actual original completion/resource
journals, finishing at 21:25 UTC. The full archive remains outside Git and its
checksum/location are recorded in the compact completion receipt. The
original design-v2 job automatically started after this closure.

## Inputs and scope

The frozen [production plan](../metadata/full_expanded_model_inputs_plan_20261002.json)
requires the closed logical case index, directed measurement catalog v3,
matched measurement join v2 and expanded covariance index. It retains the full
54-scenario matching design, both alignment guides, all four domain-matching
policies, both residue masks, six original length/coverage screens and both
mask-specific and joint-mask eligibility gates. No matching is repeated.

Ten Parquet partitions contain every case under each mask and the four physical
order combinations `00`, `01`, `10`, `11` plus their mean: **751,880 input rows**.
The first digit identifies the target alignment order and the second identifies
the control order. The mean requires both usable target orders and both usable
control orders. Eligibility continues to require the original checks across
both alignment orders even for individual order contrasts.

Every row links the original four gene endpoints, gene nodes, target/control
physical pairs, family component and species-pattern row. Shared predicted
models do not collapse distinct genes. The covariance index links the full
9,812 patterns on 526 tips and five working species-tree alternatives; older
pattern factors do not cover this expanded cohort.

The full **51,840 setting-count rows**, including empty settings, retain the
original target population and **56,965,652 unmatched decisions** through the
closed matching attrition ledger. Exclusions and same-model comparisons remain
explicit. Stable row identifiers support exact joins; they are not independent
replicates.

## Responses and predictors

Responses are target minus control CA RMSD in angstroms and target minus control
native TM dissimilarity, defined as one minus the mean of the two
endpoint-normalized native TM scores. These native endpoint scores share a
single alignment and are not independently optimized TM alignments.

Two separate sequence axes are prepared: gene-tree patristic distance and exact
identity of the structurally aligned residues. Each has first, second and third
power contrasts. For power `k`, the contrast is `target^k - control^k`.
For the mean order comparison, transform each order before averaging. Squaring
the mean identity would define a different model and is rejected by software
contracts. Identity measured in structure-derived alignments is conditional on
the structural correspondence and needs separate circularity controls.

Six nuisance contrasts capture mean endpoint coverage relative to original
full protein lengths, log aligned residue-pair count, the fraction of aligned
pairs whose two pLDDT values are at least 70, log geometric mean original
endpoint length, mean full-model CA pLDDT and the mean full-model fraction below
pLDDT 50. Confidence retains its original units: mean pLDDT on 0–100 and
fractions on 0–1. Log contrasts use natural logarithms.

Unusable, quarantined and identical-model states have null structural responses
and null alignment-derived predictors. Their independently available gene-tree
and full-protein predictors remain recorded. True numerical zero remains zero.
Missing measurements must not be interpreted as structural conservation.
All field definitions are in the [data dictionary](../metadata/full_expanded_model_inputs_data_dictionary_20261002.tsv).

## Independent verification and resources

The producer joins closed inputs with explicit Arrow types. A separate reader
reconstructs every identity, availability gate, state key and numerical cell
from raw source data using 50-digit Decimal arithmetic. Scale-aware tolerance
accounts for cancellation in powered gene-distance differences. SQLite
independently reconstructs all original selection memberships and each
setting's unique targets, controls, physical pairs, families, components,
species patterns, focal taxa and maximum reuse.

Inverse-reuse weight sums equal the number of represented groups; they are not
effective sample sizes. This handoff reports these sums but does not fit a
weighted covariance model. Original memberships remain necessary to construct
per-observation weights and dependence models.

The [software contract receipt](../metadata/full_expanded_model_input_fixture_validation_20261002.json)
records complete synthetic handoffs, interrupted replay, completed-restart
refusal and rejection of all 18 rehashed altered exports. The synthetic tests
are software checks, not a biological pilot or acceptance of production data.
Final closure requires agreement between producer and reader, all source and
artifact hashes, and the two original process completion/resource journals.

Resources were recorded before launch: **two CPU cores, 32 GiB memory, no swap,
one BLAS thread**, a 32 GiB output/scratch allowance and 100 GiB free-space
reserve. No GPU or paid resources are used. The 1–24 hour planning range per
stage is uncalibrated and is not a project ETA. Full source hashing and Decimal/
SQLite verification may dominate runtime.

```bash
/home/bizon/anaconda3/bin/python scripts/launch_full_expanded_model_inputs.py \
  --plan metadata/full_expanded_model_inputs_plan_20261002.json
```

This command already launched the frozen production version and refuses an
existing launch. Reproduction requires a new explicit output/launch identity;
do not overwrite original jobs or evidence. The
[launch inventory](../metadata/full_expanded_model_inputs_launches_20261002.json)
records producer, reader and closure handles.

## Next scientific requirements

The fixed grid implies **622,080 future model-setting records** from two
responses, two sequence axes and three polynomial degrees, or **3,110,400
nominal tree-setting fits** over five tree alternatives. These fits have not
been launched. Full-cohort design ranks, exact observation/design equivalence
and a resource estimate must precede fitting. Shared cases across settings
must not multiply effective sample size or be treated as independent tests.

Expanded covariance fitting, model adequacy, calibrated uncertainty and
multiple-testing controls remain necessary. Five working substitution kernels
do not establish accepted dated species trees or reconciled gene histories.
Prediction-source, domain, PAE, missingness, uneven sampling and ascertainment
controls remain open. The input stage alone establishes no structural
acceleration, selection, ecological association or accepted ancestral result.
