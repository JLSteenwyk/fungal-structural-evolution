#!/usr/bin/env python3
"""Run one source-verified taxon through the gated full-panel RepeatMasker stage."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def census(path: Path) -> tuple[int, int]:
    records = bases = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def write_json_atomic(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--library-builder", type=Path, required=True)
    parser.add_argument("--repeatmasker", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parent.parent
    rm_config = json.loads((root / config["repeatmodeler_config"]).read_text())
    with (root / rm_config["input_manifest"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != config["taxa"] or not 0 <= args.task_index < len(rows):
        raise ValueError("task index or complete-panel manifest is invalid")
    row = rows[args.task_index]
    source = Path(row["genome_fasta_gz"]).resolve()
    if not source.is_file() or sha256(source) != row["sha256"]:
        raise ValueError(f"source FASTA verification failed for {row['taxon_id']}")
    rm_receipt = (root / rm_config["output_root"] / row["taxon_id"] / "receipt.json").resolve()
    if not rm_receipt.is_file():
        raise FileNotFoundError(f"RepeatModeler completion unavailable for {row['taxon_id']}")

    work = (root / config["output_root"] / row["taxon_id"]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    receipt_path = work / "receipt.json"
    if receipt_path.is_file():
        existing = json.loads(receipt_path.read_text())
        outputs = [Path(p) for p in existing.get("outputs", [])]
        if (existing.get("status") == "completed" and existing.get("source_sha256") == row["sha256"]
                and outputs and all(path.is_file() and path.stat().st_size > 0 for path in outputs)):
            print(f"already completed {row['taxon_id']}")
            return
        raise RuntimeError(f"conflicting completed receipt for {row['taxon_id']}")

    combined = work / "combined-library.fasta"
    library_receipt = work / "combined-library-receipt.json"
    if not combined.exists():
        subprocess.run([
            str(args.library_builder.resolve()), "--repeatmodeler-receipt", str(rm_receipt),
            "--dfam-fasta", str((root / config["dfam_fungal_consensus"]["path"]).resolve()),
            "--dfam-sha256", config["dfam_fungal_consensus"]["sha256"],
            "--output", str(combined), "--receipt", str(library_receipt),
        ], check=True)
    library = json.loads(library_receipt.read_text())
    if sha256(combined) != library.get("combined_library_sha256"):
        raise ValueError(f"combined-library checksum mismatch for {row['taxon_id']}")

    fasta = work / "source.fasta"
    if not fasta.exists():
        with gzip.open(source, "rb") as input_handle, fasta.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
    records, bases = census(fasta)
    if records != int(row["fasta_records"]) or bases != int(row["fasta_bases"]):
        raise ValueError(f"decompressed FASTA census mismatch for {row['taxon_id']}")
    output_dir = work / "annotation"
    output_dir.mkdir(exist_ok=True)
    log = work / "run.log"
    command = [str(args.repeatmasker.resolve()), *config["repeatmasker"]["options"], "-lib", str(combined), "-dir", str(output_dir), str(fasta)]
    with log.open("a") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.flush()
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    expected = [
        output_dir / "source.fasta.out", output_dir / "source.fasta.tbl",
        output_dir / "source.fasta.gff", output_dir / "source.fasta.align",
    ]
    missing = [str(path) for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"RepeatMasker output missing or empty: {missing}")
    result = {
        "schema_version": 1,
        "status": "completed",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "taxon_id": row["taxon_id"],
        "species_name": row["species_name"],
        "assembly_accession": row["assembly_accession"],
        "source_fasta_gz": str(source),
        "source_sha256": row["sha256"],
        "source_fasta_records": records,
        "source_fasta_bases": bases,
        "combined_library": str(combined),
        "combined_library_sha256": library["combined_library_sha256"],
        "repeatmasker_version": config["repeatmasker"]["version"],
        "engine": config["repeatmasker"]["engine"],
        "command": command,
        "outputs": [str(path) for path in expected],
        "output_sha256": {str(path): sha256(path) for path in expected},
        "scope": "Per-taxon repeat-instance annotation. It does not identify rearrangements or support genome-architecture or structural-evolution inference alone."
    }
    write_json_atomic(receipt_path, result)
    print(f"completed {row['taxon_id']}")


if __name__ == "__main__":
    main()
