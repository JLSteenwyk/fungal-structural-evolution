#!/usr/bin/env python3
"""Join completed assembly and gene-order context inputs for all manifest taxa."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--resource-census", type=Path, required=True)
    parser.add_argument("--context-receipt", type=Path, required=True)
    parser.add_argument("--context-table", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.table.exists() or args.receipt.exists():
        raise FileExistsError("fresh output paths required")
    manifest = read_tsv(args.manifest)
    census = json.loads(args.resource_census.read_text())
    context_receipt = json.loads(args.context_receipt.read_text())
    context = read_tsv(args.context_table)
    if census.get("status") != "passed_receipt_level_full_panel_genome_architecture_resource_census":
        raise ValueError("completed resource census required")
    if context_receipt.get("status") != "complete_full_panel_gene_order_context_summary":
        raise ValueError("completed context summary required")
    if context_receipt.get("table_sha256") != sha256(args.context_table):
        raise ValueError("changed context table")
    taxa = {row["taxon_id"] for row in manifest}
    if len(manifest) != 526 or len(taxa) != 526 or census.get("taxa") != 526 or context_receipt.get("taxa") != 526:
        raise ValueError("full 526-taxon inputs required")
    census_rows = {row["taxon_id"]: row for row in census["taxon_dispositions"]}
    context_rows = {row["taxon_id"]: row for row in context}
    if set(census_rows) != taxa or set(context_rows) != taxa:
        raise ValueError("input taxon sets disagree")
    fields = [
        "taxon_id", "study_role", "lineage", "assembly_accession", "manifest_assembly_level",
        "busco_complete_pct", "contamination_status", "coordinate_table_status", "gene_order_rows",
        "gene_order_seqids", "singleton_seqids", "max_rows_on_seqid", "native_gene_rows",
        "provisional_orf_rows", "source_product_links", "selected_product_links", "adjacent_pairs",
        "overlapping_adjacent_pairs", "gap_pairs", "gap_bp_q25", "gap_bp_median", "gap_bp_q75",
        "gap_bp_max", "assembly_statistics_status", "receipt_assembly_level", "assembly_total_length",
        "assembly_scaffold_count", "assembly_scaffold_n50",
    ]
    args.table.parent.mkdir(parents=True, exist_ok=True)
    missing_stats = 0
    with args.table.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in manifest:
            taxon = source["taxon_id"]
            census_row, context_row = census_rows[taxon], context_rows[taxon]
            if census_row["assembly_accession"] != source["assembly_accession"]:
                raise ValueError(f"assembly mismatch: {taxon}")
            status = census_row["assembly_statistics_status"]
            missing_stats += status == "not_available"
            writer.writerow({
                "taxon_id": taxon, "study_role": source["study_role"], "lineage": source["lineage"],
                "assembly_accession": source["assembly_accession"],
                "manifest_assembly_level": source["assembly_level"], "busco_complete_pct": source["busco_complete_pct"],
                "contamination_status": source["contamination_status"],
                "coordinate_table_status": context_row["coordinate_table_status"],
                "gene_order_rows": context_row["rows"], "gene_order_seqids": context_row["seqids"],
                "singleton_seqids": context_row["singleton_seqids"], "max_rows_on_seqid": context_row["max_rows_on_seqid"],
                "native_gene_rows": context_row["native_gene_rows"], "provisional_orf_rows": context_row["provisional_orf_rows"],
                "source_product_links": context_row["source_product_links"], "selected_product_links": context_row["selected_product_links"],
                "adjacent_pairs": context_row["adjacent_pairs"], "overlapping_adjacent_pairs": context_row["overlapping_adjacent_pairs"],
                "gap_pairs": context_row["gap_pairs"], "gap_bp_q25": context_row["gap_bp_q25"],
                "gap_bp_median": context_row["gap_bp_median"], "gap_bp_q75": context_row["gap_bp_q75"], "gap_bp_max": context_row["gap_bp_max"],
                "assembly_statistics_status": status, "receipt_assembly_level": census_row.get("assembly_level", ""),
                "assembly_total_length": census_row.get("total_length", ""), "assembly_scaffold_count": census_row.get("scaffold_count", ""),
                "assembly_scaffold_n50": census_row.get("scaffold_n50", ""),
            })
    if missing_stats != census["assembly_statistics_status_counts"].get("not_available", 0):
        raise ValueError("assembly-statistics missingness count differs")
    result = {
        "schema_version": 1, "status": "complete_full_panel_genome_architecture_covariate_inputs",
        "taxa": 526, "missing_assembly_statistics_taxa": missing_stats, "table": str(args.table),
        "table_sha256": sha256(args.table), "source_hashes": {str(path): sha256(path) for path in (
            args.manifest, args.resource_census, args.context_receipt, args.context_table)},
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "A source-bound covariate-input join only. Assembly and annotation fragmentation, gaps and coordinate status remain explicit analysis covariates; this table does not call repeats, synteny, rearrangements, orthology, duplication or structural associations.",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "taxa", "missing_assembly_statistics_taxa")}, indent=2))


if __name__ == "__main__":
    main()
