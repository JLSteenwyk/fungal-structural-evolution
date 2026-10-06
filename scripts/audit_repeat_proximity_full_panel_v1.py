#!/usr/bin/env python3
"""Audit all checksum-bound repeat-instance and gene-proximity outputs."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


OUTPUT_NAMES = {
    "repeat_instances.tsv", "repeat_instances.summary.json", "gene_repeat_relations.tsv",
    "gene_repeat_summary.tsv", "gene_repeat_proximity.receipt.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def status_for(receipt_path: Path, row: dict, config: dict, root: Path) -> tuple[str, dict]:
    if not receipt_path.is_file():
        return "not_started", {}
    try:
        receipt = json.loads(receipt_path.read_text())
    except (OSError, json.JSONDecodeError):
        return "unreadable_receipt", {}
    if receipt.get("status") != "completed":
        return "noncompleted_receipt", receipt
    if receipt.get("taxon_id") != row["taxon_id"] or receipt.get("source_sha256") != row["sha256"]:
        return "manifest_identity_or_source_mismatch", receipt
    try:
        rm_receipt_path = Path(receipt["repeatmasker_receipt"])
        rm_receipt = json.loads(rm_receipt_path.read_text())
        rm_out = next(Path(value) for value in rm_receipt["outputs"] if value.endswith(".out"))
    except (KeyError, OSError, StopIteration, json.JSONDecodeError):
        return "unreadable_upstream_repeatmasker_receipt", receipt
    if (rm_receipt.get("status") != "completed" or rm_receipt.get("taxon_id") != row["taxon_id"]
            or rm_receipt.get("source_sha256") != row["sha256"]
            or receipt.get("repeatmasker_out_sha256") != rm_receipt.get("output_sha256", {}).get(str(rm_out))
            or not rm_out.is_file() or sha256(rm_out) != receipt.get("repeatmasker_out_sha256")):
        return "upstream_repeatmasker_binding_or_checksum_mismatch", receipt
    if Path(receipt.get("gene_order", "")).resolve() != (root / config["gene_order_root"] / f"{row['taxon_id']}.gene_order.tsv.gz").resolve():
        return "gene_order_path_mismatch", receipt
    output_paths = [Path(value) for value in receipt.get("outputs", [])]
    if {path.name for path in output_paths} != OUTPUT_NAMES or len(output_paths) != len(OUTPUT_NAMES):
        return "unexpected_output_set", receipt
    hashes = receipt.get("output_sha256", {})
    for path in output_paths:
        if not path.is_file() or path.stat().st_size == 0 or hashes.get(str(path)) != sha256(path):
            return "output_missing_empty_or_checksum_mismatch", receipt
    try:
        parsed = json.loads(next(path for path in output_paths if path.name == "repeat_instances.summary.json").read_text())
        linked = json.loads(next(path for path in output_paths if path.name == "gene_repeat_proximity.receipt.json").read_text())
    except (OSError, json.JSONDecodeError):
        return "unreadable_child_receipt", receipt
    if parsed.get("taxon_id") != row["taxon_id"] or linked.get("taxon_id") != row["taxon_id"]:
        return "child_receipt_taxon_mismatch", receipt
    if receipt.get("repeat_instances") != parsed.get("instances") or receipt.get("genes") != linked.get("genes"):
        return "parent_child_count_mismatch", receipt
    return "completed_verified_repeat_proximity", receipt


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
    rm_config = json.loads((root / config["repeatmasker_config"]).read_text())
    discovery = json.loads((root / rm_config["repeatmodeler_config"]).read_text())
    manifest_path = (root / discovery["input_manifest"]).resolve()
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != config["taxa"] or len({row["taxon_id"] for row in rows}) != config["taxa"]:
        raise ValueError("frozen full-panel manifest does not match configured taxon count")
    gene_readback = (root / config["gene_order_readback"]).resolve()
    readback = json.loads(gene_readback.read_text())
    if readback.get("status") != "passed_full_panel_gene_order_readback" or readback.get("taxa") != len(rows):
        raise ValueError("full gene-order readback does not establish complete-panel coverage")
    statuses, details, totals = Counter(), [], Counter()
    output_root = (root / config["output_root"]).resolve()
    for row in rows:
        status, result = status_for(output_root / row["taxon_id"] / "receipt.json", row, config, root)
        statuses[status] += 1
        if status == "completed_verified_repeat_proximity":
            totals["repeat_instances"] += result["repeat_instances"]
            totals["genes"] += result["genes"]
            totals["relations"] += sum(result.get("relation_counts", {}).values())
        details.append({"taxon_id": row["taxon_id"], "status": status})
    complete = statuses.get("completed_verified_repeat_proximity", 0) == len(rows) and len(statuses) == 1
    payload = {
        "schema_version": 1,
        "status": "passed_full_panel_repeat_proximity_audit" if complete else "full_panel_repeat_proximity_progress_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": str(config_path), "configuration_sha256": sha256(config_path),
        "manifest": str(manifest_path), "manifest_sha256": sha256(manifest_path),
        "gene_order_readback": str(gene_readback), "gene_order_readback_sha256": sha256(gene_readback),
        "taxa_expected": len(rows), "counts": dict(sorted(statuses.items())),
        "verified_totals": dict(sorted(totals.items())), "taxa": details,
        "scope": "Independent complete-panel integrity audit of output names, bytes, hashes, taxon/source identity, immutable coordinate binding, and child-summary counts. It does not establish annotation correctness, repeat burden, rearrangement, gene birth, duplication, or structural-evolution effects."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(payload["status"], payload["counts"])


if __name__ == "__main__":
    main()
