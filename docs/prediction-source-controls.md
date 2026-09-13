# Prediction-source controls and expanded mapping

Prediction pipelines are now explicit strata in marker mapping. The production default is provider **GDM**, tool **AlphaFold Monomer v2.0 pipeline**. Only exact provider/tool matches are eligible; there is no fallback to another pipeline when that source is missing. The deterministic confidence/version/model-ID ranking is applied within that source, and does not establish accuracy. The previous immutable structural snapshots remain historical checkpoints.

`map_marker_structures.py` freezes complete inventory records before mapping and audits source availability for every marker-protein link. A concurrent unfinished final record is excluded; malformed complete records raise an error. The frozen inventory is saved outside Git with a checksum. Alternatives remain in the inventory and are counted separately from source-eligible models. All 125 markers and the 526-taxon design remain in scope.

The expanded GDM mapping was launched with:

```bash
python scripts/map_marker_structures.py --output results/structural_markers/gdm-expanded-v1
```

Its terminal receipt establishes the final coverage; an existing output directory or partially written residue table does not establish completion. The previous 422-model snapshot continues to anchor the already executed comparisons until a new mapping, PAE and native-state audit are complete. An explicit provider-selection test verifies that a higher-confidence model from another tool cannot replace the requested source.

## Available same-sequence comparison

Every cross-pipeline, exact-sequence model pair in the frozen inventory was compared. At this checkpoint, there is **one distinct sequence** with such a pair: UniProt B9W6V0, 173 residues, GDM model AF-B9W6V0-F1 v6 and contributed ATBC model AF-0000000007270089 v1. The latter reports **ColabFold/AlphaFold2 Monomer v1.5.5 Pipeline**. This is the available subset of the full inventory, not a separately selected pilot or a representative sample of pipeline effects.

Full polymer identity, coordinate checksums, unambiguous complete CA positions, finite coordinates and confidence in [0,100] were verified. The literal polymer-sequence field is accepted when the canonical field is absent, with exact sequence equality and no guessing of modification codes. Both model-specific PAE matrices were retrieved and validated. The contributed PAE response is gzip-encoded; the reader now decodes transport compression before parsing and records the response checksum/encoding separately from the locally compressed JSON artifact.

At pLDDT ≥70 in both models, 171/173 residues remain. Their proper-rotation CA RMSD is **1.086 Å**, and mean PAE-filtered local CA-distance change is **0.193 Å**, retaining 3,331/3,337 eligible local pairs. At pLDDT ≥90, 141 residues remain, RMSD is 0.979 Å and local mean change is 0.163 Å. Confidence filtering changes the residue set; these differences do not show which model is more accurate.

Both pipelines use sequence-derived AlphaFold-family predictions. This comparison does not provide experimental validation, eliminate circularity, quantify systematic source effects across fungi, or establish equivalence of confidence scales. Broader cross-pipeline coverage, ESMFold controls and experimental structures remain required. In particular, do not treat pipeline-induced differences as evolutionary changes between identical sequences.

```bash
python scripts/compare_prediction_sources.py \
  --inventory results/structural_markers/gdm-expanded-v1/source_inventory.jsonl \
  --output results/structure_sources/available-pairs-v1
python -m unittest discover -s tests -p test_structure_source_selection.py
```

The comparison uses all available cross-pipeline pairs, records pLDDT thresholds 0/70/90, and separately labels insufficient coverage. Local CA pairs are within 15 Å in either model and separated by at least three sequence positions; directional PAE must be ≤10 Å in both models for the filtered metric. These are the same geometric definitions used elsewhere in the project. The one-pair resource estimate was recorded before execution: one CPU worker, cached/serial PAE retrieval, less than 1 GB RAM and 10 MB output expected, seconds to minutes, no paid resources.

Receipts and results are in `metadata/prediction_source_*`. The large frozen inventory and model/PAE data remain outside Git. A figure is not needed for this single available sequence; the complete three-row table is retained for reproducibility.
