#!/usr/bin/env python3
"""Build an independently verified full-panel genomic-DNA manifest for repeats.

The script is intentionally source-bound: it accepts no taxon selection and
requires one verified DNA receipt for every row in the supplied study manifest.
It performs a fresh SHA-256 read of every compressed FASTA before emitting a
tabular input manifest and a compact receipt.  It does not run a repeat caller.
"""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--receipt-root", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.output_manifest.exists() or args.receipt.exists():
        raise FileExistsError("refusing to overwrite an existing manifest or receipt")

    with args.study_manifest.open(newline="") as handle:
        study_rows = list(csv.DictReader(handle, delimiter="\t"))
    study_by_id = {row["taxon_id"]: row for row in study_rows}
    if len(study_by_id) != len(study_rows):
        raise ValueError("study manifest contains duplicate taxon_id values")

    receipt_paths = sorted(args.receipt_root.glob("*/receipt.json"))
    if len(receipt_paths) != len(study_rows):
        raise ValueError(
            f"expected {len(study_rows)} DNA receipts, found {len(receipt_paths)}"
        )

    out_rows = []
    total_compressed_bytes = 0
    total_bases = 0
    total_records = 0
    for receipt_path in receipt_paths:
        source = json.loads(receipt_path.read_text())
        taxon = source.get("taxon", {})
        taxon_id = taxon.get("taxon_id")
        if taxon_id not in study_by_id:
            raise ValueError(f"DNA receipt has unknown taxon_id: {taxon_id!r}")
        if source.get("status") != "verified_publisher_bound_genomic_dna":
            raise ValueError(f"unqualified DNA receipt: {taxon_id}")
        for field in ("species_name", "assembly_accession", "genome_url"):
            if taxon.get(field) != study_by_id[taxon_id].get(field):
                raise ValueError(f"manifest/receipt disagreement for {taxon_id} {field}")
        fasta = Path(source["path"])
        if not fasta.is_file():
            raise FileNotFoundError(f"DNA file unavailable for {taxon_id}: {fasta}")
        observed_sha256 = sha256(fasta)
        if observed_sha256 != source.get("sha256"):
            raise ValueError(f"DNA SHA-256 mismatch for {taxon_id}")
        stats = source.get("fasta", {})
        records = stats.get("fasta_records")
        bases = stats.get("total_length")
        if not isinstance(records, int) or not isinstance(bases, int):
            raise ValueError(f"missing FASTA census for {taxon_id}")
        total_compressed_bytes += fasta.stat().st_size
        total_bases += bases
        total_records += records
        out_rows.append({
            "taxon_id": taxon_id,
            "species_name": taxon["species_name"],
            "study_role": taxon["study_role"],
            "assembly_accession": taxon["assembly_accession"],
            "assembly_level": taxon["assembly_level"],
            "genome_fasta_gz": str(fasta),
            "sha256": observed_sha256,
            "publisher_md5": source.get("publisher_md5", ""),
            "source_url": source["url"],
            "fasta_records": records,
            "fasta_bases": bases,
            "receipt_path": str(receipt_path),
        })
    if {row["taxon_id"] for row in out_rows} != set(study_by_id):
        raise ValueError("DNA receipts do not cover precisely the study manifest taxa")

    out_rows.sort(key=lambda row: row["taxon_id"])
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.output_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)
    manifest_sha256 = sha256(args.output_manifest)
    result = {
        "schema_version": 1,
        "status": "completed_full_panel_repeat_input_readback",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "study_manifest": str(args.study_manifest),
        "study_manifest_taxa": len(study_rows),
        "receipt_root": str(args.receipt_root),
        "output_manifest": str(args.output_manifest),
        "output_manifest_sha256": manifest_sha256,
        "taxa": len(out_rows),
        "counts_by_role": {
            role: sum(row["study_role"] == role for row in out_rows)
            for role in sorted({row["study_role"] for row in out_rows})
        },
        "compressed_fasta_bytes": total_compressed_bytes,
        "fasta_records": total_records,
        "fasta_bases": total_bases,
        "scope": (
            "Fresh SHA-256 readback of every source-bound compressed genomic FASTA "
            "for the complete study panel. This is an input gate only; it does not "
            "annotate repeats or make biological inferences."
        ),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
