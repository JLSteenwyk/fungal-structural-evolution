#!/usr/bin/env python3
"""Configure the local, pinned FamDB installation to use an audited data directory."""
import argparse
import configparser
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--famdb-conf", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    if not args.famdb_conf.is_file() or not args.data_dir.is_dir():
        raise FileNotFoundError("FamDB configuration or data directory is unavailable")
    files = sorted(args.data_dir.glob("*.h5"))
    if not files:
        raise ValueError("FamDB data directory has no HDF5 partitions")
    configuration = configparser.ConfigParser()
    configuration.read(args.famdb_conf)
    if not configuration.has_section("famdb"):
        configuration.add_section("famdb")
    configuration.set("famdb", "FAMDB_DATA_DIR", str(args.data_dir.resolve()))
    with args.famdb_conf.open("w") as handle:
        configuration.write(handle)
    receipt = {
        "schema_version": 1,
        "status": "configured_local_famdb_data_directory",
        "configured_utc": datetime.now(timezone.utc).isoformat(),
        "famdb_conf": str(args.famdb_conf.resolve()),
        "data_dir": str(args.data_dir.resolve()),
        "partitions": [
            {"name": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        ],
        "scope": "Local runtime configuration for the pinned FamDB installation; it does not annotate repeats or alter source Dfam files."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
