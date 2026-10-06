#!/usr/bin/env python3
"""Export one taxon's immutable gene dispositions and selected-product anchors."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ANCHOR_FIELDS = [
    "taxon_id", "study_role", "coordinate_status", "seqid", "start", "end", "strand",
    "feature_row_id", "feature_id", "ordinal_on_seqid", "previous_gap_bp", "next_gap_bp",
    "product_id", "selected_product_count", "source_product_count",
    "product_coordinate_anchor_multiplicity", "product_coordinate_mapping_status",
]
DISPOSITION_FIELDS = [
    "taxon_id", "study_role", "coordinate_status", "seqid", "start", "end", "strand",
    "feature_row_id", "feature_id", "ordinal_on_seqid", "previous_gap_bp", "next_gap_bp",
    "source_product_count", "selected_product_count", "source_product_ids_json", "selected_product_ids_json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def read_manifest(path: Path, expected_taxa: int) -> list[dict]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != expected_taxa or len({row["taxon_id"] for row in rows}) != expected_taxa:
        raise ValueError("analysis manifest does not have the configured unique taxon count")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parent.parent
    manifest_path = (root / config["analysis_manifest"]).resolve()
    manifest = read_manifest(manifest_path, config["taxa"])
    if not 0 <= args.task_index < len(manifest):
        raise ValueError("task index is outside the frozen full-panel manifest")
    producer_path = (root / config["gene_order_receipt"]).resolve()
    producer = json.loads(producer_path.read_text())
    readback = json.loads((root / config["gene_order_readback"]).read_text())
    if (producer.get("status") != "complete_full_panel_immutable_gene_order_export"
            or producer.get("taxa") != config["taxa"]
            or readback.get("status") != "passed_full_panel_gene_order_readback"
            or readback.get("taxa") != config["taxa"]):
        raise ValueError("full immutable gene-order producer/readback gate has not passed")
    source_by_taxon = {item["taxon_id"]: item for item in producer.get("taxon_dispositions", [])}
    if set(source_by_taxon) != {row["taxon_id"] for row in manifest}:
        raise ValueError("gene-order receipt and frozen analysis manifest differ")
    manifest_row = manifest[args.task_index]
    taxon_id = manifest_row["taxon_id"]
    source = (root / source_by_taxon[taxon_id]["table"]).resolve()
    if not source.is_file() or sha256(source) != source_by_taxon[taxon_id]["table_sha256"]:
        raise ValueError("immutable gene-order table is missing or checksum-invalid")

    work = (root / config["output_root"] / taxon_id).resolve()
    work.mkdir(parents=True, exist_ok=True)
    anchors = work / "selected_product_anchors.tsv.gz"
    dispositions = work / "gene_coordinate_dispositions.tsv.gz"
    receipt_path = work / "receipt.json"
    if receipt_path.exists() or anchors.exists() or dispositions.exists():
        raise FileExistsError("refusing to overwrite gene-product anchor outputs")
    temporary_anchors = work / f".{anchors.name}.tmp"
    temporary_dispositions = work / f".{dispositions.name}.tmp"
    if temporary_anchors.exists() or temporary_dispositions.exists():
        raise FileExistsError("refusing to reuse incomplete temporary gene-product anchor output")
    counts = Counter()
    coordinate_statuses = Counter()
    product_coordinate_multiplicity = Counter()
    try:
        with gzip.open(source, "rt", newline="") as source_handle:
            reader = csv.DictReader(source_handle, delimiter="\t")
            required = set(DISPOSITION_FIELDS) - {"study_role"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError("immutable gene-order table lacks required fields")
            for row in reader:
                if row["taxon_id"] != taxon_id:
                    raise ValueError("gene-order row taxon identity differs from manifest")
                selected = json.loads(row["selected_product_ids_json"])
                source_products = json.loads(row["source_product_ids_json"])
                if not isinstance(selected, list) or not isinstance(source_products, list):
                    raise ValueError("product identifier fields must be JSON arrays")
                if len(selected) != int(row["selected_product_count"]) or len(source_products) != int(row["source_product_count"]):
                    raise ValueError("product counts differ from immutable product arrays")
                if len(set(selected)) != len(selected) or any(not isinstance(item, str) or not item for item in selected):
                    raise ValueError("selected product identifiers must be nonempty and unique per coordinate row")
                product_coordinate_multiplicity.update(selected)
        with gzip.open(source, "rt", newline="") as source_handle, gzip.open(temporary_anchors, "wt", newline="") as anchor_handle, gzip.open(temporary_dispositions, "wt", newline="") as disposition_handle:
            reader = csv.DictReader(source_handle, delimiter="\t")
            anchor_writer = csv.DictWriter(anchor_handle, fieldnames=ANCHOR_FIELDS, delimiter="\t")
            disposition_writer = csv.DictWriter(disposition_handle, fieldnames=DISPOSITION_FIELDS, delimiter="\t")
            anchor_writer.writeheader()
            disposition_writer.writeheader()
            for row in reader:
                if row["taxon_id"] != taxon_id:
                    raise ValueError("gene-order row taxon identity differs from manifest")
                selected = json.loads(row["selected_product_ids_json"])
                disposition_writer.writerow({
                    field: (manifest_row["study_role"] if field == "study_role" else row[field])
                    for field in DISPOSITION_FIELDS
                })
                counts["gene_coordinate_rows"] += 1
                coordinate_statuses[row["coordinate_status"]] += 1
                for product_id in selected:
                    multiplicity = product_coordinate_multiplicity[product_id]
                    anchor_writer.writerow({
                        "taxon_id": taxon_id, "study_role": manifest_row["study_role"],
                        "coordinate_status": row["coordinate_status"], "seqid": row["seqid"],
                        "start": row["start"], "end": row["end"], "strand": row["strand"],
                        "feature_row_id": row["feature_row_id"], "feature_id": row["feature_id"],
                        "ordinal_on_seqid": row["ordinal_on_seqid"], "previous_gap_bp": row["previous_gap_bp"],
                        "next_gap_bp": row["next_gap_bp"], "product_id": product_id,
                        "selected_product_count": row["selected_product_count"], "source_product_count": row["source_product_count"],
                        "product_coordinate_anchor_multiplicity": multiplicity,
                        "product_coordinate_mapping_status": ("unique_coordinate_anchor" if multiplicity == 1 else "multiple_coordinate_anchors"),
                    })
                    counts["selected_product_anchors"] += 1
                    counts["unique_coordinate_anchor_rows" if multiplicity == 1 else "multiple_coordinate_anchor_rows"] += 1
        os.replace(temporary_anchors, anchors)
        os.replace(temporary_dispositions, dispositions)
    except Exception:
        for temporary in (temporary_anchors, temporary_dispositions):
            if temporary.exists():
                temporary.unlink()
        raise
    receipt = {
        "schema_version": 2, "status": "completed", "completed_utc": datetime.now(timezone.utc).isoformat(),
        "taxon_id": taxon_id, "study_role": manifest_row["study_role"],
        "gene_order_table": str(source), "gene_order_table_sha256": source_by_taxon[taxon_id]["table_sha256"],
        "gene_order_producer": str(producer_path), "gene_order_readback": str((root / config["gene_order_readback"]).resolve()),
        "outputs": [str(anchors), str(dispositions)],
        "output_sha256": {str(path): sha256(path) for path in (anchors, dispositions)},
        "counts": dict(sorted(counts.items())), "coordinate_status_counts": dict(sorted(coordinate_statuses.items())),
        "scope": "All immutable coordinate rows and selected source-product links are retained. A selected product linked to multiple coordinate rows is explicitly retained as multiple_coordinate_anchors rather than arbitrarily selected or discarded. These links are annotation mappings only, not an orthology, HOG, synteny, duplication, rearrangement, or structural-evolution result."
    }
    atomic_json(receipt_path, receipt)
    print(f"completed {taxon_id}: {counts['gene_coordinate_rows']} coordinate rows; {counts['selected_product_anchors']} anchors")


if __name__ == "__main__":
    main()
