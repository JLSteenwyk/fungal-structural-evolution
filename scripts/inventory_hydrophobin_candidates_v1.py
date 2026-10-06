#!/usr/bin/env python3
"""Inventory Pfam hydrophobin candidates without inferring homology or function."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path


PROFILES = {
    "PF01185.24": {"name": "Hydrophobin", "type": "Family"},
    "PF06766.17": {"name": "Hydrophobin_2", "type": "Family"},
    "PF29785.1": {"name": "Hydrophobin_D", "type": "Domain"},
}
CHUNKS = {"PF01185.24": "chunk-036", "PF06766.17": "chunk-035", "PF29785.1": "chunk-049"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--domain-db", type=Path, required=True)
    parser.add_argument("--domain-receipt", type=Path, required=True)
    parser.add_argument("--structure-registry", type=Path, required=True)
    parser.add_argument("--registry-receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--published-receipt",
        type=Path,
        help="Optional tracked copy of the compact receipt; must not already exist.",
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    domain_receipt = json.loads(args.domain_receipt.read_text())
    registry_receipt = json.loads(args.registry_receipt.read_text())
    if not str(domain_receipt.get("status", "")).startswith("complete_"):
        raise ValueError("Domain database lacks a completed receipt")
    if not str(registry_receipt.get("status", "")).startswith("complete_"):
        raise ValueError("Structure registry lacks a completed receipt")

    catalog_rows = {row["chunk"]: row for row in read_tsv(args.catalog)}
    if len(catalog_rows) != 64:
        raise ValueError("Expected all 64 completed full-Pfam catalog shards")
    annotations = {}
    source_hashes = {
        str(args.catalog): sha256(args.catalog),
        str(args.domain_receipt): sha256(args.domain_receipt),
        str(args.registry_receipt): sha256(args.registry_receipt),
        str(args.manifest): sha256(args.manifest),
    }
    for accession, chunk in CHUNKS.items():
        row = catalog_rows.get(chunk)
        if not row:
            raise ValueError(f"Missing catalog row for {chunk}")
        path = Path(row["annotation_path"])
        if sha256(path) != row["annotation_sha256"]:
            raise ValueError(f"Changed annotation shard: {path}")
        source_hashes[str(path)] = row["annotation_sha256"]
        with gzip.open(path, "rt", newline="") as handle:
            for hit in csv.DictReader(handle, delimiter="\t"):
                if hit["pfam_accession"] == accession:
                    expected = PROFILES[accession]
                    if hit["pfam_name"] != expected["name"] or hit["pfam_type"] != expected["type"]:
                        raise ValueError(f"Unexpected metadata for {accession}")
                    if hit["hit_id"] in annotations:
                        raise ValueError("Duplicate hit identifier")
                    annotations[hit["hit_id"]] = hit

    manifest = {row["taxon_id"]: row for row in read_tsv(args.manifest)}
    if len(manifest) != 526 or set(row["study_role"] for row in manifest.values()) != {"ingroup", "outgroup"}:
        raise ValueError("Unexpected frozen analysis manifest")
    domains = sqlite3.connect(f"file:{args.domain_db}?mode=ro", uri=True)
    registry = sqlite3.connect(f"file:{args.structure_registry}?mode=ro", uri=True)
    proteins = domains.cursor()
    models = registry.cursor()
    rows = []
    for hit in annotations.values():
        sequence_id = hit["sequence_id"]
        protein_rows = proteins.execute(
            "SELECT taxon_id, protein_id, query_source FROM proteins WHERE sequence_id=? ORDER BY taxon_id, protein_id",
            (sequence_id,),
        ).fetchall()
        if not protein_rows:
            raise ValueError(f"No protein link for {sequence_id}")
        sequence_hash = sequence_id.removeprefix("S")
        model = models.execute(
            "SELECT model_id, version, length, path FROM models WHERE sequence_sha256=?", (sequence_hash,)
        ).fetchone()
        for taxon_id, protein_id, query_source in protein_rows:
            taxon = manifest.get(taxon_id)
            if not taxon:
                raise ValueError(f"Domain protein taxon absent from manifest: {taxon_id}")
            record = {
                "hit_id": hit["hit_id"],
                "pfam_accession": hit["pfam_accession"],
                "pfam_name": hit["pfam_name"],
                "pfam_clan": hit["pfam_clan"],
                "sequence_id": sequence_id,
                "protein_length": hit["protein_length"],
                "taxon_id": taxon_id,
                "species_name": taxon["species_name"],
                "study_role": taxon["study_role"],
                "lineage": taxon["lineage"],
                "protein_id": protein_id,
                "query_source": query_source,
                "alignment_start": hit["alignment_start"],
                "alignment_end": hit["alignment_end"],
                "envelope_start": hit["envelope_start"],
                "envelope_end": hit["envelope_end"],
                "hmm_coverage": hit["hmm_coverage"],
                "domain_score": hit["domain_score"],
                "domain_ga": hit["domain_ga"],
                "independent_evalue": hit["independent_evalue"],
                "posterior_accuracy": hit["posterior_accuracy"],
                "model_status": "source_model_present" if model else "no_source_model",
                "model_id": model[0] if model else "",
                "model_version": model[1] if model else "",
                "model_length": model[2] if model else "",
                "model_path": model[3] if model else "",
            }
            rows.append(record)
    if not rows:
        raise ValueError("No hydrophobin profile hits were retained")
    rows.sort(key=lambda row: (row["pfam_accession"], row["taxon_id"], row["protein_id"], row["hit_id"]))
    fields = list(rows[0])
    args.output.mkdir(parents=True)
    write_tsv(args.output / "hydrophobin_candidates.tsv", rows, fields)
    profile_counts = Counter(row["pfam_accession"] for row in rows)
    taxon_profile_counts = Counter((row["taxon_id"], row["pfam_accession"]) for row in rows)
    summary = {
        "status": "complete_provisional_hydrophobin_profile_inventory",
        "accepted_profile_hits": len(annotations),
        "taxon_protein_profile_rows": len(rows),
        "unique_sequences": len({row["sequence_id"] for row in rows}),
        "taxa": len({row["taxon_id"] for row in rows}),
        "profile_row_counts": dict(sorted(profile_counts.items())),
        "profile_taxon_counts": {
            accession: len({taxon for taxon, profile in taxon_profile_counts if profile == accession})
            for accession in PROFILES
        },
        "model_status_counts": dict(sorted(Counter(row["model_status"] for row in rows).items())),
        "source_hashes": source_hashes,
        "script_sha256": sha256(Path(__file__)),
        "artifacts": {"hydrophobin_candidates.tsv": sha256(args.output / "hydrophobin_candidates.tsv")},
        "interpretation": "All accepted matches to the three named Pfam hydrophobin profiles in the completed full-Pfam catalog are retained with protein/taxon links and source-model availability. A Pfam Family hit, including one with a source model, is not an orthogroup, a structural clade, a confirmed hydrophobin function, a gene-copy call, a structural-change estimate, or a selected evolutionary case study. Reconciled family trees, coordinate/confidence qualification, structural comparison, prediction-source controls and biological validation remain required.",
    }
    receipt_text = json.dumps(summary, indent=2) + "\n"
    (args.output / "receipt.json").write_text(receipt_text)
    if args.published_receipt:
        if args.published_receipt.exists():
            raise FileExistsError(args.published_receipt)
        args.published_receipt.parent.mkdir(parents=True, exist_ok=True)
        args.published_receipt.write_text(receipt_text)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
