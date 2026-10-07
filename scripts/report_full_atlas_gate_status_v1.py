#!/usr/bin/env python3
"""Report the immutable-receipt gates required before full-atlas clustering."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

DEFAULT_RECEIPTS = (
    Path("metadata/full_atlas_coordinate_profiles_readback_20261006_v1.json"),
    Path("metadata/full_atlas_missing_pae_20261005_v2.json"),
    Path("metadata/full_atlas_missing_pae_readback_20261006_v1.json"),
)


def read_receipt(path: Path) -> dict[str, Any]:
    item: dict[str, Any] = {"path": str(path), "exists": path.is_file()}
    if not item["exists"]:
        return item
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        item["read_error"] = f"{type(error).__name__}: {error}"
        return item
    if not isinstance(value, dict):
        item["read_error"] = "receipt JSON must be an object"
        return item
    item["status"] = value.get("status")
    item["scientific_eligibility"] = value.get("scientific_eligibility")
    return item


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--receipt", action="append", type=Path, default=[],
        help="receipt to check; may be supplied repeatedly (defaults to clustering gates)",
    )
    parser.add_argument(
        "--require-ready", action="store_true",
        help="exit nonzero unless every receipt is present, readable, and has a completed/passed status",
    )
    args = parser.parse_args()
    paths = tuple(args.receipt) or DEFAULT_RECEIPTS
    receipts = [read_receipt(path) for path in paths]
    ready = all(
        item.get("exists") and not item.get("read_error")
        and isinstance(item.get("status"), str)
        and (item["status"].startswith("completed_") or item["status"].startswith("passed_"))
        for item in receipts
    )
    print(json.dumps({"all_gates_ready": ready, "receipts": receipts}, indent=2, sort_keys=True))
    if args.require_ready and not ready:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
