#!/usr/bin/env python3
"""Independently validate the full-panel genome-architecture covariate join."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_tsv(value: object) -> str:
    return "" if value is None else str(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--resource-census", type=Path, required=True)
    parser.add_argument("--context-table", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    producer = json.loads(args.producer_receipt.read_text())
    if producer.get("status") != "complete_full_panel_genome_architecture_covariate_inputs":
        raise ValueError("completed covariate producer receipt required")
    if producer.get("table_sha256") != sha256(args.table):
        raise ValueError("covariate table hash differs")
    manifest = {row["taxon_id"]: row for row in rows(args.manifest)}
    census = {row["taxon_id"]: row for row in json.loads(args.resource_census.read_text())["taxon_dispositions"]}
    context = {row["taxon_id"]: row for row in rows(args.context_table)}
    observed = rows(args.table)
    taxa = set(manifest)
    if len(taxa) != 526 or set(census) != taxa or set(context) != taxa:
        raise ValueError("source taxon sets disagree")
    if len(observed) != 526 or {row["taxon_id"] for row in observed} != taxa:
        raise ValueError("output taxon set differs")
    missing_stats = provisional = 0
    for row in observed:
        taxon = row["taxon_id"]
        source, assembly, order = manifest[taxon], census[taxon], context[taxon]
        checks = {
            "study_role": source["study_role"], "lineage": source["lineage"],
            "assembly_accession": source["assembly_accession"], "manifest_assembly_level": source["assembly_level"],
            "busco_complete_pct": source["busco_complete_pct"], "contamination_status": source["contamination_status"],
            "coordinate_table_status": order["coordinate_table_status"], "gene_order_rows": order["rows"],
            "gene_order_seqids": order["seqids"], "singleton_seqids": order["singleton_seqids"],
            "max_rows_on_seqid": order["max_rows_on_seqid"], "native_gene_rows": order["native_gene_rows"],
            "provisional_orf_rows": order["provisional_orf_rows"], "source_product_links": order["source_product_links"],
            "selected_product_links": order["selected_product_links"], "adjacent_pairs": order["adjacent_pairs"],
            "overlapping_adjacent_pairs": order["overlapping_adjacent_pairs"], "gap_pairs": order["gap_pairs"],
            "gap_bp_q25": order["gap_bp_q25"], "gap_bp_median": order["gap_bp_median"],
            "gap_bp_q75": order["gap_bp_q75"], "gap_bp_max": order["gap_bp_max"],
            "assembly_statistics_status": assembly["assembly_statistics_status"],
            "receipt_assembly_level": as_tsv(assembly.get("assembly_level", "")),
            "assembly_total_length": as_tsv(assembly.get("total_length", "")),
            "assembly_scaffold_count": as_tsv(assembly.get("scaffold_count", "")),
            "assembly_scaffold_n50": as_tsv(assembly.get("scaffold_n50", "")),
        }
        for key, expected in checks.items():
            if row[key] != expected:
                raise ValueError(f"{key} mismatch: {taxon}")
        missing_stats += row["assembly_statistics_status"] == "not_available"
        provisional += row["coordinate_table_status"] == "provisional_orf_coordinates_not_gene_order"
    if (missing_stats, provisional) != (7, 2):
        raise ValueError("missing/provisional dispositions differ")
    result = {
        "schema_version": 1, "status": "passed_full_panel_genome_architecture_covariate_input_readback",
        "taxa": 526, "missing_assembly_statistics_taxa": missing_stats, "provisional_orf_coordinate_taxa": provisional,
        "producer_receipt_sha256": sha256(args.producer_receipt), "table_sha256": sha256(args.table),
        "source_hashes": {str(path): sha256(path) for path in (args.manifest, args.resource_census, args.context_table)},
        "script_sha256": sha256(Path(__file__)), "checked_utc": datetime.now(timezone.utc).isoformat(),
        "interpretation": "Independent exact join readback passed. It confirms covariate provenance and explicit missing/provisional dispositions, not biological genome-architecture or structure-evolution conclusions.",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "taxa", "missing_assembly_statistics_taxa", "provisional_orf_coordinate_taxa")}, indent=2))


if __name__ == "__main__":
    main()
