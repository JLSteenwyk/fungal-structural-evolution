#!/usr/bin/env python3
"""Independently audit all full-panel RepeatMasker receipts and output bytes."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_SUFFIXES = ("source.fasta.out", "source.fasta.tbl", "source.fasta.out.gff", "source.fasta.align")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def completed_status(receipt_path: Path, row: dict) -> tuple[str, dict]:
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
    library = Path(receipt.get("combined_library", ""))
    if not library.is_file() or sha256(library) != receipt.get("combined_library_sha256"):
        return "invalid_combined_library", receipt
    outputs = [Path(value) for value in receipt.get("outputs", [])]
    by_name = {path.name: path for path in outputs}
    if set(by_name) != set(REQUIRED_SUFFIXES):
        return "unexpected_output_set", receipt
    hashes = receipt.get("output_sha256", {})
    for name in REQUIRED_SUFFIXES:
        path = by_name[name]
        if not path.is_file() or path.stat().st_size == 0:
            return "missing_or_empty_output", receipt
        if hashes.get(str(path)) != sha256(path):
            return "output_checksum_mismatch", receipt
    return "completed_verified_repeatmasker_output", receipt


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
    discovery = json.loads((root / config["repeatmodeler_config"]).read_text())
    manifest_path = (root / discovery["input_manifest"]).resolve()
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != config["taxa"] or len({row["taxon_id"] for row in rows}) != config["taxa"]:
        raise ValueError("frozen full-panel manifest does not match configured taxon count")
    output_root = (root / config["output_root"]).resolve()
    statuses, details = Counter(), []
    for row in rows:
        status, result = completed_status(output_root / row["taxon_id"] / "receipt.json", row)
        statuses[status] += 1
        details.append({"taxon_id": row["taxon_id"], "status": status,
                        "repeatmasker_version": result.get("repeatmasker_version"),
                        "engine": result.get("engine")})
    complete = statuses.get("completed_verified_repeatmasker_output", 0) == len(rows) and len(statuses) == 1
    payload = {
        "schema_version": 1,
        "status": "passed_full_panel_repeatmasker_output_audit" if complete else "full_panel_repeatmasker_progress_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": str(config_path),
        "configuration_sha256": sha256(config_path),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "taxa_expected": len(rows),
        "counts": dict(sorted(statuses.items())),
        "taxa": details,
        "scope": "Independent full-panel receipt, source identity, combined-library and output-byte checksum audit. It establishes artifact completeness only; it does not validate repeat biological annotation, coverage, rearrangement, gene proximity, duplication, or structural evolution."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(payload["status"], payload["counts"])


if __name__ == "__main__":
    main()
