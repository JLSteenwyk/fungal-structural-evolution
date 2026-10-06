#!/usr/bin/env python3
"""Parse and link one fully verified RepeatMasker taxon without overwriting output."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_repeatmasker_receipt(path: Path, taxon_id: str, source_sha256: str) -> dict:
    receipt = json.loads(path.read_text())
    if receipt.get("status") != "completed" or receipt.get("taxon_id") != taxon_id:
        raise ValueError("RepeatMasker receipt is not a completed receipt for this taxon")
    if receipt.get("source_sha256") != source_sha256:
        raise ValueError("RepeatMasker receipt source checksum differs from frozen manifest")
    outputs = [Path(item) for item in receipt.get("outputs", [])]
    out_files = [item for item in outputs if item.name.endswith(".out")]
    if len(out_files) != 1:
        raise ValueError("RepeatMasker receipt must list exactly one .out output")
    output = out_files[0]
    if not output.is_file() or sha256(output) != receipt.get("output_sha256", {}).get(str(output)):
        raise ValueError("RepeatMasker .out is missing or fails its receipt checksum")
    return receipt


def checked_gene_order(path: Path, taxon_id: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"gene-order table missing: {path}")
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"taxon_id", "seqid", "start", "end", "strand", "feature_row_id", "feature_id", "coordinate_status"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("gene-order table lacks required immutable-coordinate columns")
        first = next(reader, None)
    if first is not None and first["taxon_id"] != taxon_id:
        raise ValueError("gene-order table taxon does not match frozen manifest")


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--parser-script", type=Path, required=True)
    parser.add_argument("--linker-script", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parent.parent
    rm_config = json.loads((root / config["repeatmasker_config"]).read_text())
    discovery_config = json.loads((root / rm_config["repeatmodeler_config"]).read_text())
    with (root / discovery_config["input_manifest"]).open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != config["taxa"] or not 0 <= args.task_index < len(manifest):
        raise ValueError("task index or full-panel manifest is invalid")
    gene_readback = json.loads((root / config["gene_order_readback"]).read_text())
    if gene_readback.get("status") != "passed_full_panel_gene_order_readback" or gene_readback.get("taxa") != config["taxa"]:
        raise ValueError("full gene-order readback is absent or does not cover the complete panel")
    row = manifest[args.task_index]
    taxon_id = row["taxon_id"]
    rm_receipt_path = (root / rm_config["output_root"] / taxon_id / "receipt.json").resolve()
    if not rm_receipt_path.is_file():
        raise FileNotFoundError(f"RepeatMasker completion unavailable for {taxon_id}")
    rm_receipt = checked_repeatmasker_receipt(rm_receipt_path, taxon_id, row["sha256"])
    gene_order = (root / config["gene_order_root"] / f"{taxon_id}.gene_order.tsv.gz").resolve()
    checked_gene_order(gene_order, taxon_id)

    work = (root / config["output_root"] / taxon_id).resolve()
    work.mkdir(parents=True, exist_ok=True)
    receipt_path = work / "receipt.json"
    parsed = work / "repeat_instances.tsv"
    parse_summary = work / "repeat_instances.summary.json"
    relations = work / "gene_repeat_relations.tsv"
    gene_summary = work / "gene_repeat_summary.tsv"
    proximity_receipt = work / "gene_repeat_proximity.receipt.json"
    if receipt_path.exists():
        existing = json.loads(receipt_path.read_text())
        required = [Path(item) for item in existing.get("outputs", [])]
        if (existing.get("status") == "completed" and existing.get("taxon_id") == taxon_id
                and existing.get("repeatmasker_out_sha256") == rm_receipt["output_sha256"][next(str(p) for p in rm_receipt["outputs"] if p.endswith(".out"))]
                and required and all(item.is_file() and item.stat().st_size > 0 for item in required)):
            print(f"already completed {taxon_id}")
            return
        raise RuntimeError(f"conflicting existing proximity receipt for {taxon_id}")
    if any(path.exists() for path in (parsed, parse_summary, relations, gene_summary, proximity_receipt)):
        raise RuntimeError(f"incomplete pre-existing proximity output for {taxon_id}; retain it for review and use a new stage version")

    subprocess.run([sys.executable, str(args.parser_script.resolve()), "--repeatmasker-receipt", str(rm_receipt_path),
                    "--output", str(parsed), "--summary", str(parse_summary)], check=True)
    subprocess.run([sys.executable, str(args.linker_script.resolve()), "--taxon-id", taxon_id,
                    "--gene-order", str(gene_order), "--repeat-instances", str(parsed),
                    "--relations", str(relations), "--gene-summary", str(gene_summary),
                    "--receipt", str(proximity_receipt)], check=True)
    parse_result = json.loads(parse_summary.read_text())
    link_result = json.loads(proximity_receipt.read_text())
    if parse_result.get("taxon_id") != taxon_id or link_result.get("taxon_id") != taxon_id:
        raise ValueError("downstream output taxon identity mismatch")
    out_path = next(str(path) for path in rm_receipt["outputs"] if path.endswith(".out"))
    result = {
        "schema_version": 1,
        "status": "completed",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "taxon_id": taxon_id,
        "source_sha256": row["sha256"],
        "repeatmasker_receipt": str(rm_receipt_path),
        "repeatmasker_out_sha256": rm_receipt["output_sha256"][out_path],
        "gene_order": str(gene_order),
        "gene_order_readback": str((root / config["gene_order_readback"]).resolve()),
        "outputs": [str(parsed), str(parse_summary), str(relations), str(gene_summary), str(proximity_receipt)],
        "output_sha256": {str(path): sha256(path) for path in (parsed, parse_summary, relations, gene_summary, proximity_receipt)},
        "repeat_instances": parse_result["instances"],
        "genes": link_result["genes"],
        "relation_counts": link_result["relation_counts"],
        "scope": "Checksum-bound RepeatMasker instances and immutable gene coordinates were joined without collapsing overlaps or tied nearest calls. This is coordinate evidence, not a rearrangement, causal, duplication, orthology, or structural-evolution result."
    }
    atomic_json(receipt_path, result)
    print(f"completed {taxon_id}")


if __name__ == "__main__":
    main()
