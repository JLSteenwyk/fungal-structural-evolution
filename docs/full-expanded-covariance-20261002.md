# Full expanded dependence and species covariance handoff

The full expanded cohort needs its own model inputs. The earlier 52,675-pair
covariance index and 4,568-pattern factors do not cover the expanded 75,188
logical cases. The [new producer, reader and closure](../metadata/full_expanded_covariance_launches_20261002.json)
have completed for every expanded case, including excluded and same-model cases.
The [full provenance closure](../metadata/full_expanded_covariance_completed_20261002.json)
checked 1,830,546 source/artifact bindings and both original producer/reader
completion journals. All 481,376,720 ordered covariance entries passed numerical
and original tree-edge checks. This accepts the exported working-model inputs;
biological adequacy, fitting and calibration remain separate requirements.

A complete input census finds **9,812 unique species contrast patterns** across
**3,413 families**. All endpoint taxa are present in the existing **526-tip**
working kernels. These are the expanded matching cohort, not all possible
protein comparisons across the sampling design. Every case remains linked to
its original genes, guide, node identities, all 4,250,692 selection memberships
and the complete case/mask measurement table.

## Entity reuse

The producer exports **1,052,632 endpoint/entity occurrences**: 14 per logical
case. Target and background each contribute their node, model-pair and family
identity, plus both gene endpoints and both predicted-model endpoints. A model
identity includes the model ID, integer version and coordinate checksum. Distinct
genes mapped to a shared prediction remain distinct gene entities. Repeated
endpoints of a same-model pair remain separate occurrences.

Both signed and unsigned loadings are exported. Nodes, model pairs and families
have magnitude one; gene/model endpoints have magnitude one half. Signed target
loadings are positive and background loadings negative. Unsigned loadings retain
reuse without contrast cancellation. These define prospective dependence
designs; exporting them does not select or fit a covariance model. In particular,
the signed family loading cancels when both sides belong to the same family,
whereas an unsigned family effect has a different interpretation. They must
not be substituted silently in a fitted model.

Family components connect families sharing any exported entity, across both
guides. The producer uses union-find; the separate reader constructs a sparse
family adjacency graph and computes its connected components. Both methods
recover 3,330 components across all 3,413 families. This registry
supports family blocking and reuse diagnostics. It is not a complete residual
dependence correction or proof that different components are biologically
independent. It retains repeated target nodes and physical-model reuse, which
the older background-only nesting model did not fully represent.

## Species contrast and factors

For each case the working species contrast is the focal taxon minus one half
of each background endpoint taxon. Repeated taxa are combined exactly with
integer twice-weights; shared taxa can cancel, including a completely zero
pattern. The pattern ID hashes its ordered taxon labels and twice-weights.
Its identity does not depend on the structural response, mask or matching
scenario. All 526 kernel tips remain in the column dictionary.

The complete pattern design W is decomposed on its numerical row space as
W = Q R. The rank is computed for this expanded matrix; the older rank of 242
is not imposed. The full production rank is 301. No zero pattern occurs in
this production cohort; zero-pattern behavior passed separate software
fixtures. For each working tree kernel C, a Cholesky decomposition of
R C R-transpose gives L, and F = Q L represents W C W-transpose. No jitter or
eigenvalue clipping is introduced. A failed positive-definiteness check stops
the handoff rather than inventing a factor.

The producer checks every ordered pattern-pair covariance in bounded row
blocks for each of the five current tree alternatives: **481,376,720 entries**
in total. The independent reader parses each original Newick tree with
DendroPy, constructs tip-to-edge incidence features weighted by square root
of branch length, and checks every exported covariance against the resulting
contrast path products. It also checks every source case, entity occurrence,
family component, pattern weight, factor row, rank, basis, projection and core
matrix. Zero-sum contrasts make this kernel insensitive to tree rooting; they
do not establish a root for ancestral or duplication analyses.

The kernels use substitution-distance units, not elapsed time. These five
trees remain working alternatives. Accepted species/gene phylogenies,
reconciliation, additional sensitivities and biological adequacy are separate
requirements. The factors do not estimate branch-specific structural rates or
establish that physical structural distances follow an additive tree model.

## Artifacts and reproduction

Output root: `results/phylogeny/full-expanded-covariance-20261002-v1`.
Large arrays, incidence tables and complete provenance maps stay outside Git.

| Artifact | Content |
| --- | --- |
| `case_covariance_index.tsv.gz` | One row per logical case; original identifiers, reuse count, family component, exact species-pattern row and endpoint tip indices |
| `entity_incidence.tsv.gz` | All 14 endpoint/entity occurrences per case, original entity values, signed and unsigned loadings; no endpoint collapsing |
| `family_components.tsv` | Every source family and its independently reconstructed shared-entity component |
| `patterns.tsv`, `taxa.json`, `species_contrast_design.npz` | Full taxon dictionary and exact sparse contrast design, including zero patterns |
| `pattern_basis.npy`, `species_projection.npy` | Complete numerical row-space representation |
| Five tree NPZ files | Pattern factors, reduced covariance and Cholesky factor for each original working tree |
| `receipt.json`, `readback.json`, `completion_archive.json` | Producer, separate numerical/source reader and original-journal/full-hash closure evidence |

The [field dictionary](../metadata/full_expanded_covariance_data_dictionary_20261002.tsv)
defines the tabular fields. Run the producer and reader with `--plan` pointing
to the [versioned plan](../metadata/full_expanded_covariance_plan_20261002.json).
Use a new plan/output identity to reproduce a completed run; launched plans,
scripts and accepted artifacts stay immutable. The original queued closer
writes the final completion locator after both exact original process journals
and all source/artifact hashes pass.

The [prelaunch resource estimate](../metadata/full_expanded_covariance_resources_20261002.json)
allocates two CPU equivalents, 32 GiB memory, no swap, one numerical-library
thread, 16 GiB output/scratch and a 100 GiB free-space reserve per stage.
Covariance comparisons use 128-row blocks. The uncalibrated 1–12 hour planning
allowance excludes dependencies and is not a project ETA. No GPU jobs or new
paid resources are used.

Software checks passed exact alias/reuse, same-model and zero-pattern fixtures,
interrupted full replay and completed-restart refusal. The independent reader
rejected 22 deliberately corrupted, rehashed exports. Prior case, tree and
journal proofs in these fixtures are synthetic software contracts, not
production acceptance or a biological pilot.

## Remaining model work

The expanded measurement catalog has completed independent readback and
[full provenance closure](../metadata/full_expanded_measurement_catalog_v3_completed_20261002.json)
for all 1,123,936 directed states. The complete case/mask join has passed
[separate Decimal readback and closure](../metadata/full_expanded_case_measurements_v2_completed_20261002.json)
for all 150,376 rows, with 2,369,783 source/artifact bindings and both original
journals. Complete four-order measurements exist for 70,221 full-mask and
70,015 pLDDT70 cases per outcome. These are accepted measurements, not
calibrated evolutionary effects.

Next steps are the complete expanded model-input/design census, setting-specific
entity and species factors, observation-matched model comparisons, full fitting
and calibrated uncertainty. Predictor, domain, PAE, common-residue, orientation,
ascertainment, missingness and shared-ancestry controls remain required.
Positive signed or unsigned factor loadings are model inputs, not biological
effect estimates. All eight scientific aims remain incomplete.
