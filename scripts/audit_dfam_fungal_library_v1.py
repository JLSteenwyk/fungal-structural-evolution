#!/usr/bin/env python3
"""Record a source-bound Dfam library inventory for fungal repeat analyses."""
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--famdb", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    if not args.database.is_dir() or not args.famdb.is_file():
        raise FileNotFoundError("database directory or famdb executable is unavailable")

    files = sorted(args.database.glob("*.h5"))
    if not files:
        raise ValueError("no decompressed FamDB HDF5 files found")
    info = subprocess.run(
        [str(args.famdb), "-i", str(args.database), "info"],
        check=True, text=True, capture_output=True,
    ).stdout
    fungal_check = subprocess.run(
        [str(args.famdb), "-i", str(args.database), "check", "Fungi"],
        check=True, text=True, capture_output=True,
    ).stdout
    if "Version  : 4.0" not in info or "tax id: 4751" not in fungal_check:
        raise ValueError("unexpected Dfam release or missing fungal taxonomic query")
    names = {p.name for p in files}
    required = {
        "dfam40.0.h5", "dfam40.curated.consensus.0.h5",
        "dfam40.uncurated.hmm.0.h5", "dfam40.uncurated.hmm.108.h5",
    }
    if not required.issubset(names):
        raise ValueError(f"missing required Dfam files: {sorted(required - names)}")
    receipt = {
        "schema_version": 1,
        "status": "completed_dfam_fungal_library_source_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "provider": "Dfam",
            "release": "4.0",
            "release_date": "2026-05-22",
            "url": "https://www.dfam.org/releases/current/families/FamDB/",
            "famdb_format_version": "3.0.0",
        },
        "database_directory": str(args.database),
        "files": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": digest(p)}
            for p in files
        ],
        "famdb_info": info,
        "fungal_partition_check": fungal_check,
        "scope": (
            "Records the downloaded Dfam 4.0 root and curated consensus files and "
            "the root/Fungi uncurated-HMM sensitivity partitions. It establishes source "
            "custody only; it does not annotate repeats, establish repeat identity, "
            "or support genome-architecture or structural-evolution inference."
        ),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
