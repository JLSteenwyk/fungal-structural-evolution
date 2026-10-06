#!/usr/bin/env python3
"""Audit receipt-level coordinate-anchor availability for the full panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


DIRECT_LINK = {
    "ncbi_protein_gff": "explicit_protein_id:source_product",
    "transcript_gff": "direct_transcript_parent_candidate:source_product",
    "creolimax_gtf": "explicit_gtf_transcript_candidate:source_product",
}


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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    expected = {row["taxon_id"] for row in manifest}
    if len(expected) != 526:
        raise ValueError(f"expected 526 manifest taxa, found {len(expected)}")
    completion = json.loads(args.completion.read_text())
    if completion.get("taxa") != len(expected):
        raise ValueError("completion receipt taxon count disagrees with manifest")
    receipts = {}
    for path in args.registry.glob("*.receipt.json"):
        item = json.loads(path.read_text())
        if item.get("taxon_id") in receipts:
            raise ValueError(f"duplicate receipt taxon: {item.get('taxon_id')}")
        receipts[item.get("taxon_id")] = (path, item)
    if set(receipts) != expected:
        raise ValueError("registry receipts and manifest taxa differ")

    modes, anchors, dispositions = Counter(), Counter(), []
    totals = Counter()
    for row in manifest:
        path, item = receipts[row["taxon_id"]]
        mode, features = item["mapping_mode"], item["feature_counts"]
        feature_rows, cds_rows = item["feature_rows"], features.get("CDS", 0)
        modes[mode] += 1
        totals.update(feature_rows=feature_rows, cds_rows=cds_rows, database_bytes=item["database_bytes"])
        if mode == "verified_orf_coordinates":
            if item["products_without_direct_cds_candidate"] is not None:
                raise ValueError(f"ORF mode direct-CDS mismatch: {row['taxon_id']}")
            status = "provisional_orf_coordinates_not_gene_anchors"
        else:
            if mode not in DIRECT_LINK or feature_rows <= 0 or cds_rows <= 0:
                raise ValueError(f"invalid coordinate receipt: {row['taxon_id']}")
            status = ("direct_source_product_anchor_available" if
                      item["candidate_link_counts"].get(DIRECT_LINK[mode], 0) else
                      "coordinate_records_without_direct_source_product_anchor")
        anchors[status] += 1
        dispositions.append({"taxon_id": row["taxon_id"], "study_role": row["study_role"],
                             "mapping_mode": mode, "anchor_status": status,
                             "feature_rows": feature_rows, "cds_rows": cds_rows,
                             "receipt_sha256": sha256(path)})
    result = {
        "schema_version": 1,
        "status": "passed_full_panel_genome_architecture_prerequisite_audit",
        "taxa": len(expected),
        "mapping_mode_counts": dict(sorted(modes.items())),
        "anchor_status_counts": dict(sorted(anchors.items())),
        "total_feature_rows": totals["feature_rows"],
        "total_cds_rows": totals["cds_rows"],
        "total_database_bytes": totals["database_bytes"],
        "taxon_dispositions": dispositions,
        "source_hashes": {str(args.manifest): sha256(args.manifest), str(args.completion): sha256(args.completion)},
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "Every frozen taxon has one completed registry receipt and its coordinate/anchor mode is explicit. This validates receipt-level availability, not assembly or gene-model correctness, repeat annotation, synteny blocks, orthology, rearrangement events, or structural-evolution associations.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("taxa", "mapping_mode_counts", "anchor_status_counts", "total_feature_rows", "total_cds_rows")}, indent=2))


if __name__ == "__main__":
    main()
