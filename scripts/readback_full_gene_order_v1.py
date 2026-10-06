#!/usr/bin/env python3
"""Write an independent full-panel gene-order readback receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from check_full_gene_order_v1 import validate


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    checked = validate(args.producer_receipt, args.output)
    result = {
        "schema_version": 1,
        **checked,
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "producer_receipt_sha256": sha256(args.producer_receipt),
        "reader_script_sha256": sha256(Path(__file__)),
        "checker_script_sha256": sha256(Path(__file__).with_name("check_full_gene_order_v1.py")),
        "interpretation": "Independent table/hash/schema/order/gap readback passed. This verifies export integrity, not annotation correctness, orthology, synteny, duplication, repeats or structural evolution.",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(checked, indent=2))


if __name__ == "__main__":
    main()
