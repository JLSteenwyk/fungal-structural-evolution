#!/usr/bin/env python3
"""Create a receipt-level, all-taxon resource census for genome-context work.

The census deliberately reads only compact completed registry and assembly
receipts.  It does not open the coordinate SQLite databases, infer synteny,
call repeats, or make biological eligibility claims.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--completion", type=Path, required=True)
    parser.add_argument("--assembly-statistics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gene-order-bytes-per-row", type=int, default=512)
    parser.add_argument("--context-bytes-per-row", type=int, default=384)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.gene_order_bytes_per_row <= 0 or args.context_bytes_per_row <= 0:
        raise ValueError("row-size planning assumptions must be positive")

    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    taxa = [row["taxon_id"] for row in manifest]
    if len(manifest) != 526 or len(set(taxa)) != 526:
        raise ValueError("expected exactly 526 unique manifest taxa")
    completion = json.loads(args.completion.read_text())
    if completion.get("taxa") != len(taxa):
        raise ValueError("completion receipt taxon count disagrees with manifest")

    registry = {}
    for path in args.registry.glob("*.receipt.json"):
        item = json.loads(path.read_text())
        taxon = item.get("taxon_id")
        if taxon in registry:
            raise ValueError(f"duplicate registry receipt taxon: {taxon}")
        registry[taxon] = (path, item)
    if set(registry) != set(taxa):
        raise ValueError("registry receipts and manifest taxa differ")

    assemblies = {}
    for path in args.assembly_statistics.glob("*.receipt.json"):
        item = json.loads(path.read_text())
        taxon = item.get("taxon_id")
        if taxon in assemblies:
            raise ValueError(f"duplicate assembly receipt taxon: {taxon}")
        assemblies[taxon] = (path, item)

    modes, assembly_status = Counter(), Counter()
    totals = Counter()
    dispositions = []
    for row in manifest:
        taxon = row["taxon_id"]
        receipt_path, receipt = registry[taxon]
        features = receipt["feature_counts"]
        gene_rows = features.get("gene", 0)
        selected = receipt["selected_representatives"]
        modes[receipt["mapping_mode"]] += 1
        totals.update(feature_rows=receipt["feature_rows"], cds_rows=features.get("CDS", 0),
                      gene_rows=gene_rows, selected_representatives=selected,
                      database_bytes=receipt["database_bytes"])
        assembly_item = assemblies.get(taxon)
        disposition = {
            "taxon_id": taxon,
            "study_role": row["study_role"],
            "assembly_accession": row["assembly_accession"],
            "mapping_mode": receipt["mapping_mode"],
            "feature_rows": receipt["feature_rows"],
            "gene_rows": gene_rows,
            "cds_rows": features.get("CDS", 0),
            "selected_representatives": selected,
            "registry_receipt_sha256": sha256(receipt_path),
        }
        if assembly_item is None:
            disposition["assembly_statistics_status"] = "not_available"
            assembly_status["not_available"] += 1
        else:
            assembly_path, assembly = assembly_item
            if assembly.get("assembly_accession") != row["assembly_accession"]:
                raise ValueError(f"assembly accession mismatch for {taxon}")
            if assembly.get("status") != "verified":
                raise ValueError(f"assembly receipt is not verified for {taxon}")
            metrics = assembly.get("primary_assembly_metrics", {})
            disposition.update({
                "assembly_statistics_status": "verified_and_accession_matched",
                "assembly_receipt_sha256": sha256(assembly_path),
                "assembly_level": assembly.get("headers", {}).get("Assembly level"),
                "total_length": metrics.get("total-length"),
                "scaffold_count": metrics.get("scaffold-count"),
                "scaffold_n50": metrics.get("scaffold-N50"),
            })
            assembly_status["verified_and_accession_matched"] += 1
            totals.update(assembly_bases=metrics.get("total-length", 0),
                          assembly_scaffolds=metrics.get("scaffold-count", 0))
        dispositions.append(disposition)

    representatives = totals["selected_representatives"]
    # One row per retained representative; local neighborhoods use at most a
    # left and right relation per row before any orthology/block expansion.
    planning = {
        "gene_order_rows": representatives,
        "local_directed_neighbor_rows_upper_bound": 2 * representatives,
        "gene_order_tsv_bytes_upper_bound": representatives * args.gene_order_bytes_per_row,
        "per_gene_context_tsv_bytes_upper_bound": representatives * args.context_bytes_per_row,
        "local_neighbor_tsv_bytes_upper_bound": 2 * representatives * args.context_bytes_per_row,
        "planning_assumptions": {
            "gene_order_bytes_per_row": args.gene_order_bytes_per_row,
            "context_bytes_per_row": args.context_bytes_per_row,
            "scope": "Conservative flat-text bounds for immutable coordinate/gene-context tables only; they are not estimates for repeat annotation, orthology, collinear-block candidates, alignments, structural comparisons, runtime, memory, or disk required by downstream tools.",
        },
    }
    result = {
        "schema_version": 1,
        "status": "passed_receipt_level_full_panel_genome_architecture_resource_census",
        "taxa": len(taxa),
        "assembly_statistics_status_counts": dict(sorted(assembly_status.items())),
        "mapping_mode_counts": dict(sorted(modes.items())),
        "totals": dict(sorted(totals.items())),
        "planning": planning,
        "taxon_dispositions": dispositions,
        "source_hashes": {
            str(args.manifest): sha256(args.manifest),
            str(args.completion): sha256(args.completion),
        },
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "This receipt-level census establishes the scale and compact-source availability for full-panel planning. Seven taxa lack a matched verified assembly-statistics receipt and remain explicit. It does not establish assembly correctness, repeat annotation, gene order, orthology, synteny, rearrangement, structural change, runtime or biological eligibility.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("taxa", "assembly_statistics_status_counts", "mapping_mode_counts", "totals", "planning")}, indent=2))


if __name__ == "__main__":
    main()
