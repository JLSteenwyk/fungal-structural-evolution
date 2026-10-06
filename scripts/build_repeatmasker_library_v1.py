#!/usr/bin/env python3
"""Create one verified RepeatMasker library from de-novo and Dfam consensus."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def copy_fasta(source: Path, destination, prefix: str) -> int:
    records = 0
    with source.open() as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
                destination.write(">" + prefix + line[1:])
            else:
                destination.write(line)
    if records == 0:
        raise ValueError(f"empty FASTA library: {source}")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmodeler-receipt", type=Path, required=True)
    parser.add_argument("--dfam-fasta", type=Path, required=True)
    parser.add_argument("--dfam-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.receipt.exists():
        raise FileExistsError("refusing to overwrite a combined library or receipt")
    de_novo = json.loads(args.repeatmodeler_receipt.read_text())
    library = Path(de_novo.get("classified_library", ""))
    if de_novo.get("status") != "completed" or not library.is_file():
        raise ValueError("RepeatModeler receipt is not completed with a library")
    if sha256(library) != de_novo.get("classified_library_sha256"):
        raise ValueError("RepeatModeler classified-library checksum mismatch")
    if not args.dfam_fasta.is_file() or sha256(args.dfam_fasta) != args.dfam_sha256:
        raise ValueError("Dfam consensus checksum mismatch")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    with temporary.open("w") as handle:
        de_novo_records = copy_fasta(library, handle, "RM2__")
        dfam_records = copy_fasta(args.dfam_fasta, handle, "DFAM40__")
    os.replace(temporary, args.output)
    result = {
        "schema_version": 1,
        "status": "completed_verified_repeatmasker_library",
        "taxon_id": de_novo["taxon_id"],
        "repeatmodeler_library": str(library),
        "repeatmodeler_library_sha256": de_novo["classified_library_sha256"],
        "repeatmodeler_records": de_novo_records,
        "dfam_consensus": str(args.dfam_fasta),
        "dfam_consensus_sha256": args.dfam_sha256,
        "dfam_records": dfam_records,
        "combined_library": str(args.output),
        "combined_library_sha256": sha256(args.output),
        "combined_library_bytes": args.output.stat().st_size,
        "scope": "Custom RepeatMasker search library only; no repeat-instance annotation or biological inference is made here."
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
