# Full-atlas structural clustering handoff

The current Foldseek database contains all **2,961,055** source models in the
frozen prediction atlas: 2,935,733 AFDB models and 25,322 ESMFold models. Both
predictors and every retained alternative remain indexed. These are source-model
records, not independent proteins, orthogroups, homologous families, or
evolutionary observations.

The producer and independent native-record readback are complete:

- [database producer receipt](../metadata/full_atlas_foldseek_database_20261006_v2.json)
- [independent native-record receipt](../metadata/full_atlas_foldseek_readback_20261006_v3.json)
- [database workflow and source-coordinate method](full-atlas-foldseek-database-20261006.md)

The native database is
`results/structural_clusters/full-prediction-atlas-database-20261006-v2` and
contains source aliases, amino-acid records, 3Di records, and float32 C-alpha
coordinates. The native readback checked all aliases, sequence hashes, 3Di
length/alphabet constraints, and coordinate dimensions/finiteness. It does not
independently reconstruct 3Di states or establish structural homology.

## Required gates before clustering

The current database is necessary but not sufficient for full-atlas clustering.
A full clustering launch must require all of the following immutable receipts,
then bind their hashes in its own plan:

1. The completed independent coordinate-profile readback
   (`metadata/full_atlas_coordinate_profiles_readback_20261006_v1.json`). It is
   active at this handoff and independently rechecks original coordinate,
   residue-confidence, archive, and rejection dispositions.
2. The completed full missing-AFDB PAE retrieval producer and independent
   reader (`metadata/full_atlas_missing_pae_20261005_v2.json` and
   `metadata/full_atlas_missing_pae_readback_20261006_v1.json`). Both the
   download and reader preserve unsuccessful dispositions; neither success rate
   nor an available PAE matrix is a confidence-calibration result.
3. A fresh resource preflight that measures available memory, disk, and the
   full-database temporary/output footprint immediately before native execution.
   The older 1,290,278-model clustering resource envelope must not be silently
   reused for this 2,961,055-model database.
4. A full source-identity map that preserves the taxon, representative-protein,
   predictor, and alternative-model relations for every cluster member. No
   source model may be dropped merely because it is an alternative prediction.

The completed [full identity-map receipt](../metadata/full_atlas_cluster_identity_map_20261007_v1.json)
binds all 2,961,055 database aliases to 3,020,392 preserved source relations
across all 526 taxa. It retains many-to-one mappings for identical sequence
models and the documented F87326 alternative-product ESMFold model rather than
forcing a representative-protein assignment or removing the model.

Use `python3 scripts/report_full_atlas_gate_status_v1.py` from the repository
root to inspect the three immutable-receipt gates. Its `--require-ready` mode
exits nonzero until every receipt is present, readable, and has a completed or
passed status; it is a launch guard, not evidence of biological eligibility.

After that guard passes, create a fresh host measurement immediately before the
native run; do not reuse a prior measurement:

```bash
python3 scripts/preflight_full_prediction_atlas_clustering_v1.py \
  --database results/structural_clusters/full-prediction-atlas-database-20261006-v2 \
  --output-parent results/structural_clusters \
  --receipt metadata/full_atlas_clustering_resource_preflight_YYYYMMDD_v1.json
```

The preflight checks the exact 2,961,055-row lookup, records the protected
full-database footprint, free filesystem capacity and available host memory,
and requires three times the protected database footprint as free capacity by
default. It measures resources only; the completed coordinate/PAE receipts,
launch-specific hashes and later independent cluster readback remain required.

Once all three gates and the fresh preflight pass, use
[`prepare_full_prediction_atlas_clustering_plan_v1.py`](../scripts/prepare_full_prediction_atlas_clustering_plan_v1.py)
to write the one immutable plan consumed by the clustering runner. The generator
requires the completed database producer/readback receipts, all three completed
coordinate/PAE receipts, the fresh passed preflight, the full source-identity
map and explicit Foldseek clustering arguments. It hashes every bound input and
refuses an incomplete/running receipt, a mismatched lookup count or a reused
plan/output path. The arguments remain an explicit, reviewable analytical
choice; they must not be supplied implicitly or inherited from an older atlas.

The pending MAFFT reconciliation recovery and 16-condition species-tree
sensitivity collection are separate requirements for later branch and
duplication interpretation. They are not substitutes for the coordinate/PAE
gates, and structural clusters must not be called orthogroups before
reconciliation and homology review.

## Planned native analysis and required audits

The launch will use the qualified Foldseek executable recorded by the database
workflow, with CPU-only execution and a newly pinned command. Its plan must
declare coverage, E-value, sensitivity, hit-cap, clustering mode, reassignment,
and confidence-seeding settings; retain native logs and temporary-path
provenance; and fail if any full-database lookup member is omitted or repeated.
The [gated clustering runner](../scripts/cluster_full_prediction_atlas_v1.py)
is self-tested for complete partition membership and duplicate-member rejection.
It refuses to run without launch-specific hashes for every completed prerequisite
receipt and every protected database artifact.
The separate [independent cluster reader](../scripts/readback_full_prediction_atlas_clusters_v1.py)
uses its own lookup and membership parser, rechecks the protected database after
native execution, and refuses to accept a partial partition or a producer-summary
disagreement.

The immediate result will be a *candidate structural partition*. An independent
reader must verify exact partition membership, representative self-membership,
source/taxon joins, and unchanged database artifacts after native execution.
Required follow-up analyses include threshold and predictor sensitivity,
direct/member-level alignment validation, confidence and PAE qualification,
remote-homology annotation, and family/reconciliation consistency checks.

Only after those checks can cluster-aware data contribute to the project’s
branch-specific structural-change, sequence–structure coupling, duplication,
domain-architecture, ecological-transition, or case-study analyses. Structural
similarity, inferred homology, and experimentally established function remain
distinct claims.
