#!/usr/bin/env python3
"""Audit full-panel immutable selected-product coordinate-anchor exports."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


OUTPUT_NAMES = {"selected_product_anchors.tsv.gz", "gene_coordinate_dispositions.tsv.gz"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parent.parent
    with (root / config["analysis_manifest"]).open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    producer_path = (root / config["gene_order_receipt"]).resolve()
    producer = json.loads(producer_path.read_text())
    source = {row["taxon_id"]: row for row in producer.get("taxon_dispositions", [])}
    if len(manifest) != config["taxa"] or len(source) != config["taxa"] or set(source) != {row["taxon_id"] for row in manifest}:
        raise ValueError("frozen manifest and immutable gene-order receipt do not have the configured full taxon set")
    statuses, details, totals = Counter(), [], Counter()
    for item in manifest:
        taxon = item["taxon_id"]
        receipt_path = (root / config["output_root"] / taxon / "receipt.json").resolve()
        if not receipt_path.is_file():
            status, result = "not_started", {}
        else:
            try:
                result = json.loads(receipt_path.read_text())
                outputs = [Path(value) for value in result.get("outputs", [])]
                hashes = result.get("output_sha256", {})
                expected = source[taxon]
                valid = (result.get("status") == "completed" and result.get("taxon_id") == taxon
                         and result.get("study_role") == item["study_role"]
                         and result.get("gene_order_table_sha256") == expected["table_sha256"]
                         and {path.name for path in outputs} == OUTPUT_NAMES and len(outputs) == len(OUTPUT_NAMES)
                         and all(path.is_file() and path.stat().st_size > 0 and hashes.get(str(path)) == sha256(path) for path in outputs)
                         and result.get("counts", {}).get("gene_coordinate_rows") == expected["rows"]
                         and result.get("counts", {}).get("selected_product_anchors") == expected["selected_product_links"])
                status = "completed_verified_gene_product_anchors" if valid else "invalid_receipt_or_output"
            except (OSError, json.JSONDecodeError):
                status, result = "unreadable_receipt", {}
        statuses[status] += 1
        if status == "completed_verified_gene_product_anchors":
            totals["gene_coordinate_rows"] += result["counts"]["gene_coordinate_rows"]
            totals["selected_product_anchors"] += result["counts"]["selected_product_anchors"]
        details.append({"taxon_id": taxon, "status": status})
    complete = statuses.get("completed_verified_gene_product_anchors", 0) == len(manifest) and len(statuses) == 1
    result = {
        "schema_version": 1,
        "status": "passed_full_panel_gene_product_anchor_audit" if complete else "full_panel_gene_product_anchor_progress_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": str(config_path), "configuration_sha256": sha256(config_path),
        "gene_order_receipt": str(producer_path), "gene_order_receipt_sha256": sha256(producer_path),
        "taxa_expected": len(manifest), "counts": dict(sorted(statuses.items())),
        "verified_totals": dict(sorted(totals.items())), "taxa": details,
        "scope": "Independent full-panel checks of frozen taxon identity, immutable gene-order source hashes and counts, output names, bytes and checksums. It does not validate annotation correctness or establish orthology, HOG membership, synteny, rearrangement, duplication or structural evolution."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"], result["counts"])


if __name__ == "__main__":
    main()
