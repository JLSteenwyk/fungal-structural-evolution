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
COORDINATE_STATE = Path(
    "results/structures/full-atlas-coordinate-profiles-readback-20261006-v1/state.json"
)
PAE_STATE = Path(
    "results/structures/full-atlas-missing-pae-retrieval-20261005-v2/state.json"
)
ATLAS_UNION = Path("metadata/full_prediction_atlas_union_completed_20261005_v1.json")


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


def read_json_object(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def completed_receipt_status(path: Path) -> str | None:
    """Return a completed/passed receipt status, if the immutable receipt proves one."""
    receipt = read_receipt(path)
    status = receipt.get("status")
    if (
        receipt.get("exists")
        and not receipt.get("read_error")
        and isinstance(status, str)
        and (status.startswith("completed_") or status.startswith("passed_"))
    ):
        return status
    return None


def progress_summary() -> dict[str, Any]:
    """Return advisory progress counters; these never determine gate readiness."""
    result: dict[str, Any] = {}
    union = read_json_object(ATLAS_UNION)
    coordinate = read_json_object(COORDINATE_STATE)
    pae = read_json_object(PAE_STATE)
    expected_models = None
    if union:
        counts = union.get("model_counts")
        if isinstance(counts, dict) and all(isinstance(v, int) for v in counts.values()):
            expected_models = sum(counts.values())
    if coordinate:
        counts = coordinate.get("counts")
        if isinstance(counts, dict):
            valid = sum(
                source.get("valid_models", 0)
                for source in counts.values()
                if isinstance(source, dict) and isinstance(source.get("valid_models", 0), int)
            )
            rejected = sum(
                source.get("rejected_models", 0)
                for source in counts.values()
                if isinstance(source, dict) and isinstance(source.get("rejected_models", 0), int)
            )
            coordinate_receipt = completed_receipt_status(DEFAULT_RECEIPTS[0])
            item: dict[str, Any] = {
                "stage": coordinate_receipt or coordinate.get("stage"),
                "valid_models": valid,
                "rejected_models": rejected,
            }
            if expected_models:
                item["expected_models"] = expected_models
                item["valid_fraction"] = valid / expected_models
            result["coordinate_profiles"] = item
    if pae:
        counts = pae.get("counts")
        expected = pae.get("expected_models")
        if isinstance(counts, dict) and isinstance(expected, int) and expected > 0:
            verified = counts.get("verified")
            failed = counts.get("failed")
            if isinstance(verified, int) and isinstance(failed, int):
                result["missing_pae"] = {
                    "stage": pae.get("stage"),
                    "verified": verified,
                    "failed": failed,
                    "expected_models": expected,
                    "verified_fraction": verified / expected,
                }
    return result


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
    print(json.dumps({
        "all_gates_ready": ready,
        "receipts": receipts,
        "progress_advisory_only": progress_summary(),
    }, indent=2, sort_keys=True))
    if args.require_ready and not ready:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
