#!/usr/bin/env python3
"""Summarize source-bound gene-order context for every exported taxon."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile(values: list[int], fraction: float) -> int | str:
    if not values:
        return ""
    values.sort()
    return values[round((len(values) - 1) * fraction)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--readback-receipt", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.table.exists() or args.receipt.exists():
        raise FileExistsError("fresh output paths required")
    producer = json.loads(args.producer_receipt.read_text())
    reader = json.loads(args.readback_receipt.read_text())
    if producer.get("status") != "complete_full_panel_immutable_gene_order_export":
        raise ValueError("completed producer receipt required")
    if reader.get("status") != "passed_full_panel_gene_order_readback":
        raise ValueError("completed independent readback required")
    if reader.get("producer_receipt_sha256") != sha256(args.producer_receipt):
        raise ValueError("reader is not bound to this producer receipt")
    rows = producer["taxon_dispositions"]
    if len(rows) != producer.get("taxa") or len({row["taxon_id"] for row in rows}) != producer.get("taxa"):
        raise ValueError("producer taxon dispositions disagree")

    header = ["taxon_id", "study_role", "coordinate_table_status", "rows", "seqids", "singleton_seqids",
              "max_rows_on_seqid", "native_gene_rows", "provisional_orf_rows", "source_product_links",
              "selected_product_links", "adjacent_pairs", "overlapping_adjacent_pairs", "gap_pairs",
              "gap_bp_q25", "gap_bp_median", "gap_bp_q75", "gap_bp_max"]
    args.table.parent.mkdir(parents=True, exist_ok=True)
    totals: Counter[str] = Counter()
    with args.table.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in rows:
            taxon = item["taxon_id"]
            path = args.tables / f"{taxon}.gene_order.tsv.gz"
            if not path.is_file() or sha256(path) != item["table_sha256"]:
                raise ValueError(f"table binding changed: {taxon}")
            seqids: Counter[str] = Counter(); gaps: list[int] = []
            counts: Counter[str] = Counter()
            with gzip.open(path, "rt", newline="") as source:
                table = csv.DictReader(source, delimiter="\t")
                for row in table:
                    if row["taxon_id"] != taxon:
                        raise ValueError(f"foreign row: {taxon}")
                    seqids[row["seqid"]] += 1; counts["rows"] += 1
                    counts[row["coordinate_status"]] += 1
                    counts["source_product_links"] += int(row["source_product_count"])
                    counts["selected_product_links"] += int(row["selected_product_count"])
                    if row["previous_gap_bp"] != "":
                        gap = int(row["previous_gap_bp"]); counts["adjacent_pairs"] += 1
                        if gap < 0: counts["overlapping_adjacent_pairs"] += 1
                        else: gaps.append(gap)
            record = {
                "taxon_id": taxon, "study_role": item["study_role"],
                "coordinate_table_status": item["coordinate_table_status"], "rows": counts["rows"],
                "seqids": len(seqids), "singleton_seqids": sum(v == 1 for v in seqids.values()),
                "max_rows_on_seqid": max(seqids.values(), default=0),
                "native_gene_rows": counts["native_gene_feature"],
                "provisional_orf_rows": counts["provisional_orf_coordinates_not_gene_order"],
                "source_product_links": counts["source_product_links"],
                "selected_product_links": counts["selected_product_links"],
                "adjacent_pairs": counts["adjacent_pairs"],
                "overlapping_adjacent_pairs": counts["overlapping_adjacent_pairs"],
                "gap_pairs": len(gaps), "gap_bp_q25": percentile(gaps, .25),
                "gap_bp_median": percentile(gaps, .5), "gap_bp_q75": percentile(gaps, .75),
                "gap_bp_max": percentile(gaps, 1),
            }
            writer.writerow(record)
            totals.update({key: value for key, value in counts.items() if isinstance(value, int)})
            totals.update(taxa=1, seqids=len(seqids), singleton_seqids=record["singleton_seqids"])
    result = {
        "schema_version": 1, "status": "complete_full_panel_gene_order_context_summary",
        "taxa": len(rows), "totals": dict(sorted(totals.items())), "table": str(args.table),
        "table_sha256": sha256(args.table), "producer_receipt_sha256": sha256(args.producer_receipt),
        "readback_receipt_sha256": sha256(args.readback_receipt), "script_sha256": sha256(Path(__file__)),
        "interpretation": "Source-bound coordinate order and adjacent-gap summary only. Scaffold fragmentation, overlap and gene density are retained as later covariates; they do not establish repeats, synteny, rearrangements, annotation errors, duplication, orthology or structural evolution.",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "taxa", "totals")}, indent=2))


if __name__ == "__main__":
    main()
