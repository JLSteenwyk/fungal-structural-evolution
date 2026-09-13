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
