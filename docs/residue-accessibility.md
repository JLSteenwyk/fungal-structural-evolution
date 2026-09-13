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
