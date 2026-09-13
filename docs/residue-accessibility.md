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
