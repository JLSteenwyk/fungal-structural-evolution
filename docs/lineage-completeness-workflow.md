# Complementary fungal lineage completeness assessment

The existing 526-taxon eukaryota_odb12.2 protein screen remains the common broad
baseline. A complementary full-ingroup batch is queued for all 501 fungal
entries. This is part of the full sampling workflow, not a pilot or a taxon
exclusion procedure. The 25 outgroups retain their broad eukaryotic results.

The six selected datasets are fungi, ascomycota, basidiomycota, chytridiomycota,
microsporidia and mucoromycota, each from OrthoDB 12.2 with creation date
2026-05-13. Five matching manifest groups receive their named panel; the
remaining fungal entries receive fungi_odb12.2. The choice is explicit and
reproducible, but does not establish that every sparse lineage has an ideal
reference panel. More specific comparisons can follow where warranted.

Dataset URLs and publisher MD5 values come from the frozen local BUSCO publisher
index. Archives are verified before extraction, then dataset identity, marker
counts and all extracted file hashes are recorded. The workflow uses local
paths and offline BUSCO runs, consistent with the [BUSCO user guide](https://busco.ezlab.org/busco_userguide).
A relative-path error in initial dataset receipt construction was caught;
the original downloader source is preserved, and an explicit resumed acquisition
reuses only archives whose publisher checksum matches. Final collection
verification and the QC batch are controlled by `advance_lineage_busco_qc.py`.

BUSCO 6.1.0 uses the same prepared protein inputs as the broad screen, four
threads per job and at most four simultaneous jobs. Before launch, every input
and dataset file is checked. Per-taxon receipts preserve commands, runtime,
source hashes, dataset identity and count summaries. Validation requires the
correct dataset and denominator, nonnegative categories summing to the marker
count, and agreement between complete and single/multicopy counts. Three tests
cover correct and inconsistent summaries. Completed outputs can be reused only
under the identical configuration and verified summary checksum; partial or
failed outputs require specific recovery.

The resource estimate allows 24 GB memory, 100 GB output and 4–96 hours on the
existing host. This is a planning envelope, not measured completion time. The
current stage is dataset acquisition/validation; the controller starts the
full batch only after that succeeds. Its target output is
`results/busco-lineage-v1`; no full lineage-QC completion is claimed yet.

```bash
python scripts/prepare_lineage_busco_datasets.py \
  --output data/busco-lineage-qc-v1
# Only for explicitly reviewed incomplete acquisition, with unchanged index:
python scripts/prepare_lineage_busco_datasets.py \
  --output data/busco-lineage-qc-v1 --resume
python scripts/run_lineage_busco_qc.py \
  --datasets data/busco-lineage-qc-v1 \
  --output results/busco-lineage-v1
python -m unittest discover -s tests -p test_lineage_busco_qc.py
```

Recovery percentages from different panels have different marker denominators.
Their differences cannot be read as a common-scale gain in genome quality.
Missing markers can reflect evolutionary loss/divergence, annotation problems
or assembly limitations; duplicated markers do not alone establish biological
duplication, contamination or ploidy. Subsequent interpretation must retain
lineage coverage, input conventions and these uncertainties, without uniform
automatic thresholds.

## Dataset verification completed and batch started

All six archives passed publisher checksums and extracted-file verification.
Marker counts are fungi 1,019; ascomycota 2,557; basidiomycota 1,811;
chytridiomycota 1,300; microsporidia 514; mucoromycota 1,541. The complete
file-hash manifest remains in `data/busco-lineage-qc-v1/receipt.json`, pinned by
its checksum in `metadata/lineage_busco_dataset_collection_receipt.json`.
The full 501-job batch is now running; its first jobs completed successfully.
Exact assignments and launch configuration are versioned. Final comparative
QC results remain pending.


## Full lineage QC and raw-output audit completed

All 501 fungal jobs completed successfully. The audit checked 998,331 marker
calls and 943,724 protein-hit rows against the dataset marker universes, input
protein identities, summary counts and frozen checksums. The combined table
retains all 526 taxa, including the 25 outgroups with their broad-panel scores.
No taxon was filtered. Versioned results are `metadata/busco_panel_quality.tsv`
and `metadata/lineage_busco_dataset_summary.tsv`; execution and audit receipts
are `metadata/lineage_busco_execution_receipt.json` and
`metadata/lineage_busco_audit_receipt.json`.

| Panel | Fungal entries | Markers | Median complete (%) |
|---|---:|---:|---:|
| Ascomycota | 234 | 2,557 | 95.58 |
| Basidiomycota | 157 | 1,811 | 96.69 |
| Chytridiomycota | 12 | 1,300 | 93.77 |
| Fungi fallback | 74 | 1,019 | 86.90 |
| Microsporidia | 12 | 514 | 96.30 |
| Mucoromycota | 12 | 1,541 | 95.00 |

These are within-panel summaries, not interchangeable completeness scales.
Protein-mode recovery measures annotation representation; missing and duplicated
calls still require biological and technical interpretation. Contamination,
assembly completeness and biological loss are not established by this audit.
The checker validates existing calls rather than independently repeating HMM
scoring.

```bash
python scripts/audit_lineage_busco_qc.py \
  --run results/busco-lineage-v1 \
  --datasets data/busco-lineage-qc-v1 \
  --output results/busco-lineage-audit-v1
```
