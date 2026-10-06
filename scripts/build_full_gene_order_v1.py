#!/usr/bin/env python3
"""Export immutable per-taxon gene-order tables from completed coordinate registries."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import tempfile
from collections import Counter
from pathlib import Path


HEADER = [
    "taxon_id", "coordinate_status", "seqid", "start", "end", "strand",
    "feature_row_id", "feature_id", "ordinal_on_seqid", "previous_gap_bp",
    "next_gap_bp", "source_product_count", "selected_product_count",
    "source_product_ids_json", "selected_product_ids_json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def product_map(connection: sqlite3.Connection) -> dict[str, list[tuple[str, int]]]:
    """Map native gene IDs to all source products and selected representatives."""
    mapped: dict[str, list[tuple[str, int]]] = {}
    for protein_id, mapping_json, selected in connection.execute(
        "SELECT protein_id, mapping_json, selected_representative FROM products"
    ):
        item = json.loads(mapping_json)
        gene_ids = json.loads(item.get("gene_ids_json", "[]"))
        for gene_id in gene_ids:
            mapped.setdefault(gene_id, []).append((protein_id, selected))
    return mapped


def native_rows(connection: sqlite3.Connection, taxon: str):
    genes = list(connection.execute(
        "SELECT row_id, seqid, start, end, strand, feature_id FROM features "
        "WHERE feature_type = 'gene' ORDER BY seqid, start, end, row_id"
    ))
    products = product_map(connection)
    ordinals: Counter[str] = Counter()
    for index, (row_id, seqid, start, end, strand, feature_id) in enumerate(genes):
        prior = genes[index - 1] if index and genes[index - 1][1] == seqid else None
        following = (genes[index + 1] if index + 1 < len(genes) and
                     genes[index + 1][1] == seqid else None)
        entries = sorted(products.get(feature_id or "", []))
        source = [protein for protein, _ in entries]
        selected = [protein for protein, chosen in entries if chosen]
        ordinals[seqid] += 1
        ordinal = ordinals[seqid]
        yield [
            taxon, "native_gene_feature", seqid, start, end, strand, row_id,
            feature_id or "", ordinal,
            "" if prior is None else start - prior[3] - 1,
            "" if following is None else following[2] - end - 1,
            len(source), len(selected), json.dumps(source, separators=(",", ":")),
            json.dumps(selected, separators=(",", ":")),
        ]


def provisional_rows(connection: sqlite3.Connection, taxon: str):
    rows = list(connection.execute(
        "SELECT source_ordinal, protein_id, contig, start, end, strand FROM provisional_orfs "
        "ORDER BY contig, start, end, source_ordinal"
    ))
    ordinals: Counter[str] = Counter()
    for index, (source_ordinal, protein_id, contig, start, end, strand) in enumerate(rows):
        prior = rows[index - 1] if index and rows[index - 1][2] == contig else None
        following = (rows[index + 1] if index + 1 < len(rows) and rows[index + 1][2] == contig else None)
        ordinals[contig] += 1
        ordinal = ordinals[contig]
        yield [
            taxon, "provisional_orf_coordinates_not_gene_order", contig, start, end, strand,
            source_ordinal, protein_id, ordinal,
            "" if prior is None else start - prior[4] - 1,
            "" if following is None else following[3] - end - 1,
            1, 0, json.dumps([protein_id]), "[]",
        ]


def write_table(connection: sqlite3.Connection, taxon: str, mode: str, destination: Path) -> dict[str, int]:
    rows = provisional_rows(connection, taxon) if mode == "verified_orf_coordinates" else native_rows(connection, taxon)
    counts: Counter[str] = Counter()
    with destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with __import__("io").TextIOWrapper(zipped, newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(HEADER)
                for record in rows:
                    writer.writerow(record)
                    counts["rows"] += 1
                    counts["source_product_links"] += int(record[11])
                    counts["selected_product_links"] += int(record[12])
                    counts["seqids"] += record[8] == 1
    return dict(counts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--completion", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--expected-taxa", type=int, default=526)
    args = parser.parse_args()
    if args.output_directory.exists() or args.receipt.exists():
        raise FileExistsError("output directory or receipt already exists")
    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != args.expected_taxa or len({x["taxon_id"] for x in manifest}) != args.expected_taxa:
        raise ValueError("manifest count does not match expected taxa")
    completion = json.loads(args.completion.read_text())
    if completion.get("taxa") != args.expected_taxa:
        raise ValueError("completion receipt does not match expected taxa")

    receipts = {}
    for path in args.registry.glob("*.receipt.json"):
        item = json.loads(path.read_text())
        if item.get("taxon_id") in receipts:
            raise ValueError(f"duplicate registry receipt: {item.get('taxon_id')}")
        receipts[item["taxon_id"]] = (path, item)
    expected = {x["taxon_id"] for x in manifest}
    if set(receipts) != expected:
        raise ValueError("registry and manifest taxa differ")

    args.output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=args.output_directory.name + ".", dir=args.output_directory.parent))
    try:
        totals: Counter[str] = Counter()
        modes: Counter[str] = Counter()
        dispositions = []
        for row in manifest:
            taxon = row["taxon_id"]
            receipt_path, item = receipts[taxon]
            database = Path(item["database"])
            if not database.is_file():
                raise FileNotFoundError(database)
            mode = item["mapping_mode"]
            destination = temporary / f"{taxon}.gene_order.tsv.gz"
            connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
            try:
                counts = write_table(connection, taxon, mode, destination)
            finally:
                connection.close()
            if mode == "verified_orf_coordinates":
                expected_rows = item["source_products"]
            else:
                expected_rows = item["feature_counts"].get("gene", 0)
            if counts.get("rows", 0) != expected_rows:
                raise ValueError(f"row count mismatch for {taxon}: {counts.get('rows', 0)} != {expected_rows}")
            modes[mode] += 1
            totals.update(counts)
            dispositions.append({
                "taxon_id": taxon, "study_role": row["study_role"], "mapping_mode": mode,
                "coordinate_table_status": ("provisional_orf_coordinates_not_gene_order" if mode == "verified_orf_coordinates"
                                            else "native_gene_features_exported"),
                "rows": counts.get("rows", 0), "seqids": counts.get("seqids", 0),
                "source_product_links": counts.get("source_product_links", 0),
                "selected_product_links": counts.get("selected_product_links", 0),
                "registry_receipt_sha256": sha256(receipt_path),
                "table": str(args.output_directory / destination.name), "table_sha256": sha256(destination),
            })
        os.replace(temporary, args.output_directory)
        temporary = None
        result = {
            "schema_version": 1,
            "status": "complete_full_panel_immutable_gene_order_export",
            "taxa": args.expected_taxa, "mapping_mode_counts": dict(sorted(modes.items())),
            "totals": dict(sorted(totals.items())), "taxon_dispositions": dispositions,
            "source_hashes": {str(args.manifest): sha256(args.manifest), str(args.completion): sha256(args.completion)},
            "script_sha256": sha256(Path(__file__)),
            "interpretation": "Every manifest taxon has a source-bound coordinate-order table. Native rows describe annotation gene features; the two ORF-coordinate taxa remain explicitly provisional. Product links are source mappings, not orthology, synteny, gene-copy or evolutionary-event calls.",
        }
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({k: result[k] for k in ("taxa", "mapping_mode_counts", "totals")}, indent=2))
    finally:
        if temporary is not None and temporary.exists():
            shutil.rmtree(temporary)


if __name__ == "__main__":
    main()
