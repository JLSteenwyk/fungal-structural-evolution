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

## Execution and remaining integration

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
with an adjacent `-audit` directory after completion. Running calculations are
not completed exposure measurements. The audit verifies source identity,
residue/atom accounting, confidence and totals; it does not independently
recompute the geometric ASA algorithm.

After completion, the reusable and new outputs still need an exact full-cohort
union, projection to the recovered paired masks, normalization, and independent
readback. The exposure/tree integration must resolve each marker through the
recovered fit collection's source inventory. Older 124-marker exposure tables
cannot substitute for the current 125-marker collection. Core/surface coupling
and functional interpretations remain downstream; extant accessibility is not
an ancestral estimate, pocket annotation or experimentally established function.
