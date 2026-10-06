#!/usr/bin/env python3
"""Independently validate a full-panel immutable gene-order export."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

from build_full_gene_order_v1 import HEADER


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    if receipt.get("status") != "complete_full_panel_immutable_gene_order_export":
        raise ValueError("unexpected receipt status")
    rows = receipt.get("taxon_dispositions", [])
    if len(rows) != receipt.get("taxa") or len({item["taxon_id"] for item in rows}) != receipt.get("taxa"):
        raise ValueError("receipt taxon disposition mismatch")
    totals: Counter[str] = Counter()
    for item in rows:
        path = args.output / f'{item["taxon_id"]}.gene_order.tsv.gz'
        if not path.is_file() or sha256(path) != item["table_sha256"]:
            raise ValueError(f"table hash mismatch: {item['taxon_id']}")
        prior_by_seqid = {}
        previous_record = None
        ordinal_by_seqid: Counter[str] = Counter()
        counts: Counter[str] = Counter()
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != HEADER:
                raise ValueError(f"table schema mismatch: {item['taxon_id']}")
            for record in reader:
                if record["taxon_id"] != item["taxon_id"]:
                    raise ValueError(f"taxon mismatch: {item['taxon_id']}")
                seqid, start, end = record["seqid"], int(record["start"]), int(record["end"])
                if start > end:
                    raise ValueError(f"inverted coordinate: {item['taxon_id']}")
                ordinal_by_seqid[seqid] += 1
                if int(record["ordinal_on_seqid"]) != ordinal_by_seqid[seqid]:
                    raise ValueError(f"ordinal mismatch: {item['taxon_id']}")
                prior = prior_by_seqid.get(seqid)
                expected_gap = "" if prior is None else str(start - prior - 1)
                if record["previous_gap_bp"] != expected_gap:
                    raise ValueError(f"previous gap mismatch: {item['taxon_id']}")
                if previous_record is not None:
                    if previous_record["seqid"] == seqid:
                        next_gap = str(start - int(previous_record["end"]) - 1)
                    else:
                        next_gap = ""
                    if previous_record["next_gap_bp"] != next_gap:
                        raise ValueError(f"next gap mismatch: {item['taxon_id']}")
                prior_by_seqid[seqid] = end
                previous_record = record
                counts.update(rows=1, source_product_links=int(record["source_product_count"]),
                              selected_product_links=int(record["selected_product_count"]))
        if previous_record is not None and previous_record["next_gap_bp"] != "":
            raise ValueError(f"terminal next gap mismatch: {item['taxon_id']}")
        counts["seqids"] = len(ordinal_by_seqid)
        for key in ("rows", "seqids", "source_product_links", "selected_product_links"):
            if counts[key] != item[key]:
                raise ValueError(f"{key} mismatch: {item['taxon_id']}")
        totals.update(counts)
    if dict(sorted(totals.items())) != receipt["totals"]:
        raise ValueError("receipt totals mismatch")
    print(json.dumps({"status": "passed_full_panel_gene_order_readback", "taxa": receipt["taxa"],
                      "totals": dict(sorted(totals.items()))}, indent=2))


if __name__ == "__main__":
    main()
