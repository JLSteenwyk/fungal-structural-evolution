#!/usr/bin/env python3
"""Link one taxon's repeat instances to immutable gene coordinates.

This produces all overlapping gene--repeat pairs and all ties for nearest
upstream/downstream calls. It never selects an arbitrary repeat among ties.
"""
import argparse
import bisect
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


RELATION_FIELDS = [
    "taxon_id", "feature_row_id", "feature_id", "coordinate_status", "seqid",
    "gene_start", "gene_end", "gene_strand", "relation", "distance_bp",
    "repeat_id", "repeat_name", "repeat_class_family", "library_origin",
    "repeat_start", "repeat_end", "repeat_strand", "repeat_overlap_marker",
]
GENE_FIELDS = [
    "taxon_id", "feature_row_id", "feature_id", "coordinate_status", "seqid",
    "gene_start", "gene_end", "gene_strand", "overlap_repeat_count",
    "upstream_nearest_distance_bp", "upstream_nearest_tie_count",
    "downstream_nearest_distance_bp", "downstream_nearest_tie_count",
    "repeat_coordinate_evidence_status",
]


def int_or_blank(value):
    return "" if value is None else int(value)


def parse_instances(path: Path, taxon_id: str):
    by_seqid = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["taxon_id"] != taxon_id:
                raise ValueError("repeat table taxon_id does not match requested taxon")
            row["query_start"] = int(row["query_start"])
            row["query_end"] = int(row["query_end"])
            by_seqid[row["query_sequence"]].append(row)
    for rows in by_seqid.values():
        rows.sort(key=lambda row: (row["query_start"], row["query_end"], int(row["repeat_id"])))
    return by_seqid


def relation(gene, repeat, kind, distance):
    return {
        "taxon_id": gene["taxon_id"],
        "feature_row_id": gene["feature_row_id"],
        "feature_id": gene["feature_id"],
        "coordinate_status": gene["coordinate_status"],
        "seqid": gene["seqid"],
        "gene_start": gene["start"], "gene_end": gene["end"], "gene_strand": gene["strand"],
        "relation": kind, "distance_bp": distance,
        "repeat_id": repeat["repeat_id"], "repeat_name": repeat["repeat_name"],
        "repeat_class_family": repeat["repeat_class_family"], "library_origin": repeat["library_origin"],
        "repeat_start": repeat["query_start"], "repeat_end": repeat["query_end"],
        "repeat_strand": repeat["strand"], "repeat_overlap_marker": repeat["overlap_marker"],
    }


def read_genes(path: Path, taxon_id: str):
    by_seqid = defaultdict(list)
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["taxon_id"] != taxon_id:
                raise ValueError("gene-order table taxon_id does not match requested taxon")
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            by_seqid[row["seqid"]].append(row)
    for rows in by_seqid.values():
        rows.sort(key=lambda row: (row["start"], row["end"], int(row["feature_row_id"])))
    return by_seqid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxon-id", required=True)
    parser.add_argument("--gene-order", type=Path, required=True)
    parser.add_argument("--repeat-instances", type=Path, required=True)
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--gene-summary", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if any(path.exists() for path in (args.relations, args.gene_summary, args.receipt)):
        raise FileExistsError("refusing to overwrite proximity outputs or receipt")
    genes = read_genes(args.gene_order, args.taxon_id)
    repeats = parse_instances(args.repeat_instances, args.taxon_id)
    args.relations.parent.mkdir(parents=True, exist_ok=True)
    args.gene_summary.parent.mkdir(parents=True, exist_ok=True)
    relation_counts = Counter()
    total_genes = total_repeats = 0
    with args.relations.open("w", newline="") as relation_handle, args.gene_summary.open("w", newline="") as gene_handle:
        relation_writer = csv.DictWriter(relation_handle, fieldnames=RELATION_FIELDS, delimiter="\t")
        gene_writer = csv.DictWriter(gene_handle, fieldnames=GENE_FIELDS, delimiter="\t")
        relation_writer.writeheader()
        gene_writer.writeheader()
        for seqid, gene_rows in genes.items():
            total_genes += len(gene_rows)
            repeat_rows = repeats.get(seqid, [])
            total_repeats += len(repeat_rows)
            starts = [row["query_start"] for row in repeat_rows]
            by_end = sorted(repeat_rows, key=lambda row: (row["query_end"], row["query_start"], int(row["repeat_id"])))
            ends = [row["query_end"] for row in by_end]
            active, repeat_cursor = [], 0
            for gene in gene_rows:
                while repeat_cursor < len(repeat_rows) and repeat_rows[repeat_cursor]["query_start"] <= gene["end"]:
                    active.append(repeat_rows[repeat_cursor])
                    repeat_cursor += 1
                active = [row for row in active if row["query_end"] >= gene["start"]]
                overlaps = sorted(active, key=lambda row: (row["query_start"], row["query_end"], int(row["repeat_id"])))
                for repeat in overlaps:
                    relation_writer.writerow(relation(gene, repeat, "overlap", 0))
                    relation_counts["overlap"] += 1
                before = bisect.bisect_left(ends, gene["start"])
                upstream = []
                upstream_distance = None
                if before:
                    nearest_end = ends[before - 1]
                    upstream_distance = gene["start"] - nearest_end
                    index = before - 1
                    while index >= 0 and ends[index] == nearest_end:
                        upstream.append(by_end[index])
                        index -= 1
                    for repeat in upstream:
                        relation_writer.writerow(relation(gene, repeat, "upstream_nearest", upstream_distance))
                        relation_counts["upstream_nearest"] += 1
                after = bisect.bisect_right(starts, gene["end"])
                downstream = []
                downstream_distance = None
                if after < len(repeat_rows):
                    nearest_start = starts[after]
                    downstream_distance = nearest_start - gene["end"]
                    index = after
                    while index < len(repeat_rows) and starts[index] == nearest_start:
                        downstream.append(repeat_rows[index])
                        index += 1
                    for repeat in downstream:
                        relation_writer.writerow(relation(gene, repeat, "downstream_nearest", downstream_distance))
                        relation_counts["downstream_nearest"] += 1
                evidence = "repeat_calls_present_on_seqid" if repeat_rows else "no_repeat_calls_on_seqid_or_seqid_unresolved"
                gene_writer.writerow({
                    "taxon_id": gene["taxon_id"], "feature_row_id": gene["feature_row_id"],
                    "feature_id": gene["feature_id"], "coordinate_status": gene["coordinate_status"],
                    "seqid": seqid, "gene_start": gene["start"], "gene_end": gene["end"],
                    "gene_strand": gene["strand"], "overlap_repeat_count": len(overlaps),
                    "upstream_nearest_distance_bp": int_or_blank(upstream_distance),
                    "upstream_nearest_tie_count": len(upstream),
                    "downstream_nearest_distance_bp": int_or_blank(downstream_distance),
                    "downstream_nearest_tie_count": len(downstream),
                    "repeat_coordinate_evidence_status": evidence,
                })
    receipt = {
        "schema_version": 1,
        "status": "completed_gene_repeat_proximity_link",
        "taxon_id": args.taxon_id,
        "gene_order": str(args.gene_order), "repeat_instances": str(args.repeat_instances),
        "relations": str(args.relations), "gene_summary": str(args.gene_summary),
        "genes": total_genes, "repeat_calls_on_gene_seqids": total_repeats,
        "relation_counts": dict(sorted(relation_counts.items())),
        "scope": "Coordinate relations only. It retains ties and overlap calls, and does not infer repeat causation, gene birth, rearrangement, orthology, or structural change."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
