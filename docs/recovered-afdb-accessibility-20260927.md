# Recovered AlphaFold accessibility — September 27, 2026

The recovered marker collection contains 30,618 models. Exact inventory checking
identified 13,127 models with reusable audited solvent-accessibility output and
17,491 models requiring new calculations (9,061,916 residues). Another 26 models
belong only to the earlier cohort and remain explicitly excluded from this union.

For every reused model, all provenance fields match the previous model except
an added `source_record_id` equal to its UniProt accession. Coordinate and sequence
hashes, version, length and paths are unchanged. Reused entry receipts and residue
table bytes were checked against the full earlier audit. An independent inventory
readback checked the exact partition, all model fields, disposition lengths and
reused entry hashes. This inventory check does not independently recalculate ASA.

`scripts/prepare_recovered_accessibility_inventory.py --plan
metadata/recovered_afdb_accessibility_inventory_plan_20260927.json` reproduces
the inventory in a fresh output directory. The output
`results/structural_annotations/accessibility-gap-afdb-recovered-20260927-v1`
contains the new-model provenance manifest, complete disposition ledger and
hash-bound receipt. It is an accessibility input subset, not a paired marker map.
See [inventory verification](../metadata/recovered_afdb_accessibility_inventory_completed_20260927.json).

## Calculation and audit (completed)

New calculations use the existing `annotate_predicted_accessibility.py` with
960 sphere points per atom and a 1.4 Å probe, matching the previous cohort.
All heavy atoms of the isolated predicted chain contribute to occlusion;
low-confidence atoms are retained and confidence is reported. Full-chain
orientation and prediction uncertainty can affect these measurements.

The pinned controller `scripts/advance_recovered_accessibility.py --config
metadata/recovered_afdb_accessibility_plan_20260927.json` runs calculation then
the full raw-coordinate output audit. It uses four CPU workers, 16 GiB RAM,
no swap and a 16 GiB output allowance, with 32 GiB available-memory and 100 GiB
free-disk preflight gates. A systematic sample of 132 prior model runtimes gives
an approximate eight-hour calculation extrapolation at four workers. The broader
4–48 hour planning allowance accommodates size effects, contention and auditing;
neither is a guaranteed ETA. GPU prediction remains paused.

Launch identity and enforced resource limits are recorded in
`metadata/recovered_afdb_accessibility_launch_20260927.json`. New ASA output is
`results/structural_annotations/accessibility-afdb-recovered-gap-20260927-v1`,
with an adjacent `-audit` directory. Both the calculation and raw-source audit
completed successfully on September 27. The audit verifies source identity,
residue/atom accounting, confidence and totals; it does not independently
recompute the geometric ASA algorithm.

The reusable and new outputs now have an exact full-cohort union, projection
to the recovered paired masks, normalization, and independent readback.
Exposure/tree integration resolves each marker through the recovered fit
collection's source inventory. Older 124-marker exposure tables
cannot substitute for the current 125-marker collection. Core/surface coupling
and functional interpretations remain downstream; extant accessibility is not
an ancestral estimate, pocket annotation or experimentally established function.

## Completed projection and normalization

The separate pinned controller waited on the calculation/audit controller's
PID, creation time and exact command and completed these stages after the
producer and full raw-source audit receipts passed:

1. Merge exactly 30,618 selected models (15,960,692 residues), using unchanged
   source entries through relative symlinks and explicitly excluding the 26
   earlier-only models.
2. Project accessibility onto all observed cells of the 125-marker paired input
   cohort: **9,453,757 residue observations**, counted directly from the paired
   amino-acid FASTAs. Masked cells remain missing.
3. Independently check every projected row against raw accessibility tables,
   mapping records and both paired alphabets.
4. Normalize using both existing reference scales, retain terminal exclusions,
   and independently check all row values and normalization summaries.

`merge_recovered_accessibility.py` accepts only the documented additional source
identifier; every other shared model field must agree. It checks all selected
coordinate, entry and residue-table hashes. CLI fixtures passed valid subset
selection and rejected missing/duplicate cohorts, undeclared exclusions, altered
residue tables, changed source identifiers/sequences and incomplete audits.
These fixtures test provenance and file handling, not ASA geometry.

The completed plan is
`metadata/recovered_afdb_accessibility_projection_plan_20260927.json`; the cohort
manifest, fixture receipt and process launch record are adjacent. The controller
uses one CPU, 48 GiB RAM, no swap and an 8 GiB output allowance, with 64 GiB
available-memory and 100 GiB free-disk gates. The 0.5–24 hour downstream planning
range is broad and does not include waiting for accessibility calculations.
The allowance accounts for the 11,202,520-row source mapping retained in memory.

The union, paired projection and normalized tables use the
`accessibility-union-afdb-recovered-20260927-v1`,
`paired-accessibility-afdb-recovered-20260927-v1`, and
`normalized-accessibility-afdb-recovered-20260927-v1` directories under
`results/structural_annotations/`. Full readbacks are written beneath
`results/recovery-20260927/afdb-accessibility-projection/`.
All five integration stages and their declared artifacts passed the completion
provenance check. Sequence–structure coupling remains downstream.

## Verified gene-tree sources and completed site summaries

The recovered collection combines 95 unchanged native marker fits with 30
refitted marker fits. `prepare_collection_topologies.py` resolved all 125 through
that collection's audited source inventory and exported relative symlinks to the
exact native AA trees. Every tree was independently checked against its native
receipt, selected source path, expanded alignment hashes and exact tip set.
The full grid contains 47,529 sites and 9,453,757 observed residue cells.
Evidence is `metadata/recovered_afdb_collection_topologies_completed_20260927.json`.
This exports existing topologies; it does not infer new trees or uncertainty.

Reproduce with `scripts/prepare_collection_topologies.py` and arguments:

- `--inputs results/phylogeny/paired-inputs-afdb-recovered-20260925-v1`
- `--models data/structural_models/garg-hochberg-v3`
- `--collection results/phylogeny/paired-fit-tables-afdb-recovered-20260926-v1`
- `--inventory metadata/recovered_afdb_fit_source_inventory_20260926.json`
- `--output results/phylogeny/collection-topologies-afdb-recovered-20260927-v1`

A separate controller completed after accessibility projection, normalization
and both readbacks. It applied the established minimum-change recurrence to AA
and 3Di states on each verified topology and summarized extant accessibility
quantiles. The independent set-based recurrence checked all 95,058 character
scores; a separate quantile implementation checked all 285,174 exposure
quantiles, including unavailable values. All observed identities, states and
model counts are checked as well.

Scripts `summarize_collection_site_exposure.py` and
`readback_collection_site_exposure.py` use the explicit topology collection;
the numerical algorithms match the existing single-run implementations. The
plan `metadata/recovered_afdb_site_exposure_plan_20260927.json` pins scripts,
inputs and the completed topology readback. The launch record identifies the
successfully completed controller. Resources are one CPU, 48 GiB RAM, no swap and 2 GiB
output, with 64 GiB available-memory and 100 GiB disk gates; the 0.25–12 hour
planning allowance excludes predecessor wait time.

These summaries and their full independent readback are complete. Parsimony
counts are alphabet-dependent minimum
changes, not rates or branch assignments. Exposure quantiles summarize extant
observations, not ancestral exposure or a phylogenetically adjusted effect.
The subsequent rate/exposure join and controlled coupling inference remain open.

## September 27 completion evidence

The calculation, projection/normalization and site-summary services are all
terminal with exit status zero. The completion verifier checked nine stages,
74 unique pinned/source/artifact files, all declared stage counts, and the
controller prerequisite hashes. Evidence is
[the completion record](../metadata/recovered_afdb_accessibility_completed_20260927.json).
Recheck from the repository root with a fresh output path:

```bash
python scripts/verify_recovered_accessibility_completion.py --output /tmp/recovered-accessibility-completion.json
```

The union has **30,618 models and 15,960,692 residues**. The paired projection
actually uses **30,342 models**, supplying **9,453,757 observations** across
**125 markers and 47,529 alignment sites**. The difference in model counts
reflects the paired projection, not additional failed accessibility calculations.
Every projected row passed raw-source readback; every normalized row passed
value and summary readback. There were no terminal-residue normalization exclusions.

Normalization retains values above one without clipping. The theoretical Tien
scale produces three such observations; the Miller scale produces 27,402.
At a relative-accessibility threshold of 0.25, 5,226,483 observations fall below
threshold using Tien and 4,827,457 using Miller; **399,026 observations (4.22%)**
change classification between reference scales. These are repeated taxon/site
observations, not independent evolutionary replicates or biological truth labels.
The downstream analysis should retain both normalization scales and continuous
accessibility, alongside confidence and domain-orientation sensitivities.

This verifies source accounting and numerical readbacks, not an independent
recalculation of the ASA geometry or equivalence to DSSP. Isolated predicted
chains omit partners; predicted domain placement can change accessibility.
Extant accessibility does not establish ancestral burial, interfaces or pockets.
The rate/exposure frame controller remains gated on the full rate optimization
and its audit; no acceleration, coupling significance or selection claim follows
from this milestone.
