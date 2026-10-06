#!/usr/bin/env python3
"""Parse one checksum-verified RepeatMasker .out into an instance table."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


FIELDS = [
    "taxon_id", "query_sequence", "query_start", "query_end", "query_left",
    "strand", "repeat_name", "repeat_class_family", "repeat_start", "repeat_end",
    "repeat_left", "repeat_id", "overlap_marker", "library_origin", "sw_score",
    "perc_divergence", "perc_deletion", "perc_insertion",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parenthesized_int(value: str) -> int:
    if not (value.startswith("(") and value.endswith(")")):
        raise ValueError(f"expected parenthesized coordinate: {value}")
    return int(value[1:-1])


def library_origin(repeat_name: str) -> str:
    if repeat_name.startswith("RM2__"):
        return "RepeatModeler_2.0.9"
    if repeat_name.startswith("DFAM40__"):
        return "Dfam_4.0_curated_consensus"
    if repeat_name.startswith("(") and repeat_name.endswith(")n"):
        return "RepeatMasker_simple_repeat"
    return "unattributed_repeatmasker_label"


def parse_row(line: str, taxon_id: str) -> dict:
    parts = line.split()
    if len(parts) < 15:
        raise ValueError(f"unexpected RepeatMasker .out row: {line.rstrip()}")
    overlap = parts[-1] == "*"
    core = parts[:-1] if overlap else parts
    if len(core) != 15:
        raise ValueError(f"unexpected RepeatMasker .out field count: {line.rstrip()}")
    score, divergence, deletion, insertion = core[:4]
    query, begin, end, left, strand, repeat, classification, repeat_begin, repeat_end, repeat_left, identifier = core[4:]
    return {
        "taxon_id": taxon_id,
        "query_sequence": query,
        "query_start": int(begin),
        "query_end": int(end),
        "query_left": parenthesized_int(left),
        "strand": strand,
        "repeat_name": repeat,
        "repeat_class_family": classification,
        "repeat_start": int(repeat_begin),
        "repeat_end": int(repeat_end),
        "repeat_left": parenthesized_int(repeat_left),
        "repeat_id": int(identifier),
        "overlap_marker": overlap,
        "library_origin": library_origin(repeat),
        "sw_score": int(score),
        "perc_divergence": float(divergence),
        "perc_deletion": float(deletion),
        "perc_insertion": float(insertion),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmasker-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.summary.exists():
        raise FileExistsError("refusing to overwrite parsed output or summary")
    receipt = json.loads(args.repeatmasker_receipt.read_text())
    if receipt.get("status") != "completed":
        raise ValueError("RepeatMasker receipt is not completed")
    output_hashes = receipt.get("output_sha256", {})
    out_files = [Path(path) for path in receipt.get("outputs", []) if str(path).endswith(".out")]
    if len(out_files) != 1 or not out_files[0].is_file():
        raise ValueError("receipt lacks exactly one accessible .out file")
    out_path = out_files[0]
    if sha256(out_path) != output_hashes.get(str(out_path)):
        raise ValueError("RepeatMasker .out checksum mismatch")
    rows = []
    with out_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("SW", "score", "There")):
                continue
            if not line[0].isdigit() and not line.lstrip()[0].isdigit():
                continue
            rows.append(parse_row(line, receipt["taxon_id"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    class_counts = Counter(row["repeat_class_family"] for row in rows)
    origin_counts = Counter(row["library_origin"] for row in rows)
    summary = {
        "schema_version": 1,
        "status": "completed_repeatmasker_instance_parse",
        "taxon_id": receipt["taxon_id"],
        "repeatmasker_receipt": str(args.repeatmasker_receipt),
        "source_out": str(out_path),
        "source_out_sha256": output_hashes[str(out_path)],
        "instance_table": str(args.output),
        "instance_table_sha256": sha256(args.output),
        "instances": len(rows),
        "classes": dict(sorted(class_counts.items())),
        "library_origins": dict(sorted(origin_counts.items())),
        "overlap_marked_instances": sum(row["overlap_marker"] for row in rows),
        "scope": "Lossless per-instance RepeatMasker parsing; overlapping calls remain explicit and no coverage union, repeat burden, rearrangement, or biological association is inferred."
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
