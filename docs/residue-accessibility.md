# Predicted residue solvent accessibility

This stage supports the planned localization of evolutionary changes in protein
cores and surfaces. It computes continuous solvent-accessible surface area on
full isolated predicted chains, preserving residue identity and CA pLDDT.
It does not assign catalytic sites, pockets, interaction interfaces or confirmed
buried/exposed biological states.

The first production run covers the frozen expanded AlphaFold snapshot:
13,153 models and 6,910,765 residues. It uses Biopython 1.86 Shrake–Rupley,
a 1.4 Å solvent probe and 960 sphere points per atom. Atomic radii, package
versions, implementation source hash, prediction snapshot and script hashes are
pinned in the run configuration. See the [Biopython SASA documentation](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html)
for the algorithm and parameters; the run uses its pinned installed version,
not whichever version is current at that URL.

Only a single predicted protein chain/model with standard C/N/O/S heavy atoms
and unambiguous residue identities is accepted. Sequence hashes and contiguous
residue numbering must match model provenance. Every residue retains its area,
CA pLDDT and atom count. Low-confidence atoms remain present during the surface
calculation because deleting them would alter occlusion. Downstream confidence
and domain-orientation sensitivity remain necessary. Missing biological
partners and ligands mean this is isolated-chain accessibility, not accessibility
in a biological complex. No residue-size normalization or categorical core
threshold is imposed at this stage.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/annotate_predicted_accessibility.py --snapshot results/structural_markers/gdm-expanded-v1 --output results/structural_annotations/accessibility-gdm-v1 --workers 4 --points 960
```

The run is already active. Per-model compressed tables and atomic completion
receipts support restart with the same pinned configuration. Do not launch a
second copy or change the live producer. Large tables remain outside Git.

Planning reserved four CPU workers with one thread each, 32 GB memory and 10 GB
output on the existing authorized host, with a deliberately broad 4–96 hour
forecast pending production timings. No new paid resources were provisioned.

Validation includes an analytic isolated-carbon sphere test, reduction in total
exposure relative to two independent spheres after adding an occluding atom,
and residue/atom area sum agreement. The first eight completed model tables
also passed checksum, full residue-grid, canonical sequence and total-area
readback. This does not independently validate the entire surface calculation
or establish numerical convergence for every protein. Sampling-resolution
sensitivity and downstream evolutionary integration remain pending.

## ESMFold coverage and numerical-resolution sensitivity

The same unmodified producer now runs on the frozen local ESMFold snapshot:
5,121 models, 1,183,341 residues, two CPU workers, 16 GB memory and 5 GB output
planning. The broad forecast is 2–24 hours. All local input proteins are within
the current 512-residue prediction limit. The first eight completed local tables
passed independent hash, residue-grid, sequence and total-area readback. Local
and AlphaFold outputs remain source-specific, and no cross-source averaging or
new branch estimates are produced by this annotation stage.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/annotate_predicted_accessibility.py --snapshot results/structural_markers/esmfold-partial-v1 --output results/structural_annotations/accessibility-esmfold-v1 --workers 2 --points 960
```

Separate resolution assessments are now running on both snapshots. Selection
uses two smallest salted model-ID hashes within fixed length bins: 1–128,
129–256, 257–512, 513–1024 and >1024 residues. AlphaFold contributes ten models;
ESMFold contributes six because its snapshot lacks proteins above 512 residues.
Selection precedes and does not depend on measured surface area. Every chosen
model is calculated at both 960 and 3,840 sphere points per atom, preserving
paired residue identities and reporting median, 95th-percentile and maximum
absolute per-residue area differences and whole-chain totals. The finer grid
is a numerical comparator, not exact surface-area ground truth.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/assess_accessibility_resolution.py --snapshot results/structural_markers/gdm-expanded-v1 --output results/structural_annotations/accessibility-resolution-gdm-v1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/assess_accessibility_resolution.py --snapshot results/structural_markers/esmfold-partial-v1 --output results/structural_annotations/accessibility-resolution-esmfold-v1
```

These assessments use one CPU each, 16 GB total memory and 2 GB total output
planning, with a 0.1–8 hour forecast reflecting length variation. They run on the
existing host without new paid resources. Both are already active; do not
launch duplicates. Their results will quantify sampling-resolution sensitivity
on these selected proteins, not prove convergence or biological validity across
the full atlas. Production annotation and evolutionary integration remain
incomplete.


## Completed local resolution assessment

The six predetermined ESMFold models (1,234 residues) completed at both sphere
resolutions. All selected-model identities, pinned files, per-residue sequence
and index grids, area totals and reported numerical summaries passed an
independent readback. Across these residues, absolute differences between 960
and 3,840 points have median 0.320 Å², 95th percentile 1.137 Å² and maximum
2.431 Å². These are residue-weighted descriptive differences in six selected
proteins, not bounds on the entire local dataset or on structural prediction
error.

Three residues switch between exactly zero and nonzero area across the two
resolutions. This supports retaining continuous area and avoiding an unqualified
zero-area definition of a biological core. The higher-resolution calculation is
still approximate. The AlphaFold assessment, including longer proteins, remains
running; its results must be considered separately.

```bash
python scripts/audit_accessibility_resolution.py --assessment results/structural_annotations/accessibility-resolution-esmfold-v1 --snapshot results/structural_markers/esmfold-partial-v1 --output results/structural_annotations/accessibility-resolution-esmfold-audit-v1
```


## Completed AlphaFold resolution assessment

All ten predetermined AlphaFold models completed at both resolutions, including
the final 1,468-residue protein. All 5,839 residue pairs passed selection,
configuration, receipt/hash, sequence/grid, total-area and numerical-summary
readback. Comparing 960 versus 3,840 sphere points gives a residue-weighted
median absolute difference of 0.282 Å², 95th percentile 1.048 Å² and maximum
2.049 Å². Thirty-nine residues change between exactly zero and nonzero area.

| Source | Selected models | Residues | Median absolute difference (Å²) | 95th percentile (Å²) | Zero/nonzero changes |
| --- | ---: | ---: | ---: | ---: | ---: |
| AlphaFold | 10 | 5,839 | 0.282 | 1.048 | 39 |
| ESMFold | 6 | 1,234 | 0.320 | 1.137 | 3 |

These source subsets differ in length and protein identity. The table is a
numerical-resolution assessment, not a comparison of predictor accuracy.
Neither set establishes convergence of every production structure. Continuous
ASA remains the primary annotation; any later normalized exposure or categorical
core/surface partition requires explicit conventions and sensitivity analysis.
Biological partners, uncertain domain orientations and low-confidence structure
remain separate sources of uncertainty.

```bash
python scripts/audit_accessibility_resolution.py --assessment results/structural_annotations/accessibility-resolution-gdm-v1 --snapshot results/structural_markers/gdm-expanded-v1 --output results/structural_annotations/accessibility-resolution-gdm-audit-v1
```


## Production residue-table audit

`scripts/audit_predicted_accessibility.py` validates production annotations
against the original mmCIF atom tables using `MMCIF2Dict`, separately from the
producer's structure-object parser. It checks source/configuration/output hashes,
sequence identity, the complete residue grid, atom names/counts, CA confidence,
finite nonnegative ASA and reported totals. This does not independently
recalculate ASA or establish biological exposure.

The default requires every model and a matching full production receipt. For
ongoing production, `--allow-partial` freezes the set of atomic entry receipts
present at audit start. Its output identifies exactly which entries were checked
and cannot establish completion of the full dataset. Unfinished entries are
left to the producer; the auditor does not modify its outputs.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/audit_predicted_accessibility.py --assessment results/structural_annotations/accessibility-esmfold-v1 --snapshot results/structural_markers/esmfold-partial-v1 --output results/structural_annotations/accessibility-esmfold-entry-audit-v1 --allow-partial
```

The first production entry audit passed all 4,652 frozen completed models and
1,059,944 residues. The other 469 models were not in this audit snapshot; this
is not full production completion. The checked model/receipt list is recorded
in `metadata/esmfold_accessibility_audited_models.tsv`, with aggregate scope and
provenance in `metadata/esmfold_accessibility_entry_audit_receipt.json`.


## Automatic full-snapshot audit

The full local-source audit is queued with
`scripts/advance_accessibility_audit.py`. Its recorded producer PID and Linux
start ticks distinguish the current calculation from any later PID reuse.
The controller waits for that producer to exit, then requires its completion
receipt and unchanged pinned configuration/code before launching the auditor
without partial mode. A stopped producer without a receipt fails explicitly.
The auditor then checks every entry and full completion counts; the controller
records success only for a full-snapshot audit result.

The resource estimate is one CPU, 4 GB RAM and 1 GB disk headroom, with 3–15
minutes expected after production completes, based on the completed entry audit.
The config is `metadata/esmfold_accessibility_audit_controller_config.json`.
This is a queued validation stage, not evidence that full validation has passed.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/advance_accessibility_audit.py --config metadata/esmfold_accessibility_audit_controller_config.json
```


## Full local production and evolutionary-site projection completed

All 5,121 frozen ESMFold models completed: 1,183,341 residues. The automatic
controller then completed the full raw-coordinate audit, covering every model
and residue. The production, audit and controller receipts are recorded in
`metadata/esmfold_accessibility_full_receipt.json`,
`metadata/esmfold_accessibility_full_audit_receipt.json` and
`metadata/esmfold_accessibility_audit_controller_receipt.json`.
This supersedes the earlier partial-audit coverage limit for this snapshot.
ASA values themselves were not independently recalculated by this audit.

`scripts/link_paired_sites_accessibility.py` now joins continuous focal-residue
ASA to every observed site in the exact paired AA/3Di inputs: 714,936 sites,
72 markers, 4,669 marker–taxon combinations and 4,630 distinct models. The
observed-site count matches the prior feature-overlap audit. Masked alignment
cells are omitted from the residue table and counted in each marker summary.
Species-specific observations are retained even when coordinates are identical.

The projection requires a complete matching accessibility audit. A rejection
check confirmed that the earlier partial audit cannot create projection outputs.
After production, all 714,936 rows passed separate readback for complete observed
site-set equality, duplicate exclusion, AA/3Di states, matrix columns and exact
ASA/confidence/atom-count/context fields. All marker cell totals were reconciled.
The original residue mapping and numerical ASA calculation are separate checks.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/link_paired_sites_accessibility.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --snapshot results/structural_markers/esmfold-partial-v1 --accessibility results/structural_annotations/accessibility-esmfold-v1 --audit results/structural_annotations/accessibility-esmfold-full-audit-v1 --output results/structural_annotations/paired-accessibility-esmfold-v1
```

The compressed residue table is outside Git under the output directory. Its
checksum, source lineage and all-marker summary are recorded in
`metadata/esmfold_paired_accessibility_receipt.json`,
`metadata/esmfold_paired_accessibility_marker_summary.tsv` and
`metadata/esmfold_paired_accessibility_readback.json`. Resources were estimated
before launch in `metadata/esmfold_paired_accessibility_resource_plan.json`.

This supplies a coordinate system for structural localization analyses. ASA
is not normalized across amino-acid types; no categorical core/surface,
interface or evolutionary transition is assigned here. Full-chain domain
placement, low-confidence occluding residues and missing biological partners
remain relevant. Phylogenetic modeling of changes in exposure is still pending.


## Reference normalization and threshold sensitivity

Normalized all 714,936 paired sites with two amino-acid-specific reference
scales: Tien et al. theoretical ALLOWED-region values and the Miller values
reproduced in their Table 1. The [primary paper](https://doi.org/10.1371/journal.pone.0080635)
recommends the theoretical scale and separates terminal residues. Its maxima
were evaluated with DSSP; our ShrakeRupley values use the same 1.4 Å probe but
algorithm/radius equivalence has not been established. These outputs are
reference-normalized accessibility indices, not a calibration of our ASA method.

`config/accessibility_normalization.json` contains all 40 denominators, source
XML URL/checksum and policies. The values were extracted from publisher XML and
all matched Biopython's Wilke and Miller tables. The normalization retains
original fields and values above one; terminal rows would receive blank indices.
No terminal residues occur among these confidence-qualified native-context sites.

| Strict threshold | Sites classified differently between scales | Fraction |
| --- | ---: | ---: |
| <0.05 | 8,300 | 1.16% |
| <0.10 | 13,147 | 1.84% |
| <0.20 | 23,845 | 3.34% |
| <0.25 | 31,332 | 4.38% |
| <0.50 | 63,609 | 8.90% |

These thresholds are diagnostic choices, not accepted biological labels.
Tien normalization produces no values above one in this site set; Miller
produces 2,846, which are retained. Absence of overshoots does not establish
method equivalence or accuracy. Changing denominators deterministically changes
threshold assignments; these counts are not independent evolutionary events.
Continuous values remain primary, with amino-acid identity and prediction
context requiring controls in later comparative models.

All original fields and 714,936 normalized rows passed readback: multiplying
each index by its scale denominator reproduced ASA within 1e-12 relative and
absolute tolerance. Both scale summaries and all five threshold summaries were
independently accumulated and matched. Receipts, summaries and validation scope
are in `metadata/esmfold_accessibility_normalization_*` and
`metadata/esmfold_accessibility_threshold_sensitivity.tsv`.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/normalize_paired_accessibility.py --projection results/structural_annotations/paired-accessibility-esmfold-v1 --snapshot results/structural_markers/esmfold-partial-v1 --config config/accessibility_normalization.json --output results/structural_annotations/paired-accessibility-normalized-esmfold-v1
```


## Topology-aware site variation and extant exposure

`scripts/summarize_site_parsimony_exposure.py` joins all 16,571 paired columns
across 72 local markers to unit-cost minimum AA and 3Di changes on their fitted
AA gene-tree topologies. Unknown tips permit any state. A postorder dynamic
program supports hard multifurcations without inventing resolutions. This is a
fixed-tree parsimony diagnostic ([method reference](https://doi.org/10.1137/0128004));
it ignores branch lengths, time, unequal substitution costs and topology uncertainty.

Each site retains observed taxon and distinct-model counts, matrix coordinates,
copy-review status, distinct character states, and median/interquartile exposure
under both normalization scales. All 714,936 taxon-site observations are accounted
for. Exposure summaries concern extant observed taxa, not reconstructed ancestors.

There are 12,759 AA-variable and 9,377 3Di-variable columns. Among the latter,
1,595 have invariant focal AA characters. Context residues, spatial partners,
other sequence differences and prediction variation can alter a 3Di state;
this does not establish structural evolution independent of sequence.
Alphabet-specific minimum counts are lower bounds, not comparable physical
units or branch-rate estimates. No branch assignments or coupling significance
are inferred. Exposure/rate models still require phylogenetic, sampling,
confidence and predictor controls.

Validation covered 243 exhaustive two-state/missing tip patterns on three tiny
topologies, plus hard-polytomy, multicolumn and invalid-input tests. For the
actual data, a separate set-membership frequency recurrence reproduced all
33,142 AA/3Di scores. Every site observation/model count and six exposure
quantiles passed readback from normalized residue rows. This validates the
calculation, not the biological adequacy of unit-cost parsimony.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/summarize_site_parsimony_exposure.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --audit results/phylogeny/paired-fit-audit-esmfold-v1 --exposure results/structural_annotations/paired-accessibility-normalized-esmfold-v1 --projection results/structural_annotations/paired-accessibility-esmfold-v1 --review config/marker_orthology_review.json --output results/phylogeny/site-parsimony-exposure-esmfold-v1
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_site_parsimony.py
```

The complete site table remains outside Git. Its hash, topology pins, source
lineage, validation and descriptive counts are recorded in
`metadata/esmfold_site_parsimony_exposure_{receipt,readback,summary}.json`;
resource planning is in `metadata/site_parsimony_exposure_resource_plan.json`.


## Full bootstrap-topology sensitivity running

`scripts/assess_site_parsimony_topologies.py` evaluates both original paired
alignments on all 1,000 saved AA UFBoot trees per marker: 72,000 topologies over
all 72 markers. It preserves repeated bootstrap trees and checks each tree's
taxon universe. Character columns remain unchanged, isolating topology choice
from the separate paired column-resampling run. Source tree archives, input
receipts, code and the parsimony helper are pinned.

Per-marker restartable outputs retain both complete bootstrap-by-site score
arrays and per-site ranges, quantiles, distinct-score counts and the fraction
matching the original-tree score. Each completed array is read back before its
marker completion receipt is written. The full stage requires every marker.
Quantiles are descriptive topology sensitivity conditional on the AA bootstrap
procedure, not calibrated confidence intervals or posterior probabilities.
No branch assignments or ancestral sequences are inferred.

Eight single-thread workers run on the existing authorized host, with 16 GB RAM
and 5 GB disk planning and a 0.25–8 hour runtime forecast. This is a launched
stage, not yet a completed all-marker result. Resource and execution configs are
`metadata/site_parsimony_topology_resource_plan.json` and
`metadata/esmfold_site_parsimony_topology_config.json`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/assess_site_parsimony_topologies.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --diagnostic results/phylogeny/site-parsimony-exposure-esmfold-v1 --output results/phylogeny/site-parsimony-topologies-esmfold-v1 --workers 8 --bootstrap-trees 1000
```


## Completed topology sensitivity and full artifact audit

All 72,000 saved bootstrap topologies completed over 16,571 sites. The retained
arrays contain 33,142,000 AA/3Di scores. Minimum counts vary across these trees
at 7,514 AA sites and 4,965 3Di sites. This records topology sensitivity rather
than a calibrated uncertainty interval or evidence of accelerated evolution.
Different alphabets and amounts of observed variation prevent interpreting
these totals as comparative accuracy or physical stability.

`scripts/audit_site_parsimony_topologies.py` checked every array's dimensions,
integer bounds and hash, all 72,000 tree tip universes, preservation of original
site annotations, every summary statistic and aggregate count. Quantiles were
reconstructed from sorted values with explicit interpolation. The audit also
selected one bootstrap tree per marker by an identity-based SHA256 rule and
independently recomputed all 33,142 site/alphabet scores on those 72 trees using
a separate set-membership recurrence. Every check passed. The remaining 999
trees per marker were not independently rescored; numerical and provenance
checks cover their stored arrays and summaries.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/audit_site_parsimony_topologies.py --assessment results/phylogeny/site-parsimony-topologies-esmfold-v1 --diagnostic results/phylogeny/site-parsimony-exposure-esmfold-v1 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fits results/phylogeny/paired-marker-fits-esmfold-partial-v1 --output results/phylogeny/site-parsimony-topologies-audit-esmfold-v1
```

Complete arrays and site tables remain outside Git. Summary/provenance and
validation artifacts are `metadata/esmfold_site_parsimony_topology_receipt.json`,
`metadata/esmfold_site_parsimony_topology_audit_receipt.json` and
`metadata/esmfold_site_parsimony_topology_marker_audit.tsv`.
The separate column-resampling run is still active. Neither analysis includes
all sources of model, alignment, prediction or biological uncertainty.

The next model-based layer is now available: [conditional site-specific evolutionary rates](site-specific-evolutionary-rates.md), covering all 72 markers and four model specifications. Controlled exposure/rate analyses remain pending.


## Full expanded AlphaFold accessibility calculation completed

The full calculation finished for all 13,153 models and 6,910,765 residues.
Its producer exited zero and the final receipt enumerates all 13,153 per-model
receipts. A compact versioned record preserves the complete receipt's path,
checksum and counts in `metadata/gdm_full_accessibility_execution_receipt.json`;
the full per-model receipt map remains outside Git with the data.

The full raw-output audit is now running. It checks all source coordinate and
output hashes, complete residue identities, heavy-atom counts, CA confidence,
finite nonnegative ASA values and per-model ASA totals against raw mmCIF tables.
This is the established audit used for the local snapshot; it does not
independently recompute solvent-accessible area or establish biological exposure,
interfaces or categorical core/surface labels. The calculations describe isolated
predicted chains, so partner binding and domain-orientation uncertainty remain
relevant.

Resource allowance: one CPU, 8 GB memory, 1 GB output and 0.2–8 planning hours on
the existing host. Plan: `metadata/gdm_full_accessibility_audit_resource_plan.json`.
Audit target: `results/structural_annotations/accessibility-gdm-full-audit-v1`.
Projection onto exact paired sites and normalization remain subsequent steps.

These two steps are now queued behind the exact running full-audit process.
`scripts/advance_accessibility_projection.py` checks the completed audit's
13,153-model, 6,910,765-residue universe and matching source receipts before
launching the existing projection and normalization programs sequentially.
The configuration pins 16 script, receipt and normalization-source files:
`metadata/gdm_accessibility_projection_controller_config.json`.
The projection covers all 124 ready markers with their existing paired masks.
Normalization uses the established Tien and Miller references without clipping;
reference thresholds remain diagnostic conventions, not validated biological
core/surface assignments.

Planning allowances are one CPU, 16 GiB memory, 10 GiB output and 0.1–8 hours;
the controller requires 32 GiB available memory before launch. These are
estimates, allowing for the 4,895,692-row mapping retained by the projection.
Outputs will be `paired-accessibility-gdm-expanded-v1` and
`paired-accessibility-normalized-gdm-expanded-v1` under
`results/structural_annotations/`. Independent projection/normalization row
readback and phylogenetically controlled exposure analyses remain pending.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/advance_accessibility_projection.py \
  --config metadata/gdm_accessibility_projection_controller_config.json
```

```bash
OPENBLAS_NUM_THREADS=1 python scripts/audit_predicted_accessibility.py \
  --assessment results/structural_annotations/accessibility-gdm-v1 \
  --snapshot results/structural_markers/gdm-expanded-v1 \
  --output results/structural_annotations/accessibility-gdm-full-audit-v1
```


## Follow-on ESMFold accessibility running

The complete follow-on cohort now has an accessibility calculation running for
4,252 models and 1,342,046 residues. All model IDs are disjoint from the previous
5,121-model accessibility cohort. The same implementation, 960 sphere points
per atom and 1.4-A probe are used; all heavy atoms remain potential occluders,
regardless of confidence. Paired-site interpretation will use the separately
qualified confidence masks.

The earlier cohort consumed 15,842.76 summed worker seconds for 1,183,341 residues.
Linear scaling by residue count gives a 1.25-hour four-worker scenario, not a
runtime bound. The resource plan allows 0.5–12 hours, four CPUs, 8 GiB memory and
5 GiB output on the existing host. See
`metadata/esmfold_followon_accessibility_resource_plan.json` and the actual
configuration in `metadata/esmfold_followon_accessibility_launch.json`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/annotate_predicted_accessibility.py \
  --snapshot results/structural_markers/esmfold-followon-complete-v1 \
  --output results/structural_annotations/accessibility-esmfold-followon-v1 \
  --workers 4 --points 960
```

The established full-output auditor is queued after the exact producer process,
using `metadata/esmfold_followon_accessibility_audit_controller_config.json`.
It pins 11 dependencies and will check raw coordinate/residue identity, hashes,
atom counts, confidence and ASA totals. This does not independently integrate
ASA. Projection, normalization and controlled exposure/evolution tests remain
pending; calculation launch is not evidence of a biological association.
