#!/usr/bin/env python3
"""Audit full-panel RepeatModeler receipts without altering any worker output."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def classified_fasta_status(path: Path) -> tuple[int, int]:
    """Return total FASTA records and records bearing a RepeatMasker class."""
    records = classified = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
                if "#" in line[1:].split(maxsplit=1)[0]:
                    classified += 1
    return records, classified


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
    with (root / config["input_manifest"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != config["taxa"]:
        raise ValueError("manifest does not match configured full-panel taxon count")
    output_root = (root / config["output_root"]).resolve()
    details = []
    statuses = Counter()
    for row in rows:
        work = output_root / row["taxon_id"]
        receipt_path = work / "receipt.json"
        if receipt_path.is_file():
            result = json.loads(receipt_path.read_text())
            library = Path(result.get("classified_library", ""))
            valid_receipt = (result.get("status") == "completed" and result.get("taxon_id") == row["taxon_id"]
                             and result.get("source_sha256") == row["sha256"] and library.is_file()
                             and sha256(library) == result.get("classified_library_sha256"))
            records, classified = classified_fasta_status(library) if valid_receipt else (0, 0)
            if valid_receipt and records > 0 and records == classified:
                status = "completed_verified_classified_library"
            elif valid_receipt:
                status = "invalid_classified_library_format"
            else:
                status = "invalid_completed_receipt"
        elif work.is_dir():
            status = "incomplete_worker_output"
        else:
            status = "not_started"
        statuses[status] += 1
        details.append({"taxon_id": row["taxon_id"], "status": status,
                        "classified_library_records": records if receipt_path.is_file() else 0,
                        "class_bearing_records": classified if receipt_path.is_file() else 0})
    payload = {
        "schema_version": 1,
        "status": "full_panel_repeatmodeler_progress_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": str(config_path),
        "taxa_expected": len(rows),
        "counts": dict(sorted(statuses.items())),
        "taxa": details,
        "scope": (
            "Independent receipt/hash audit for the full 526-taxon RepeatModeler "
            "stage. It reports execution disposition only and makes no repeat or "
            "genome-architecture inference."
        ),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
