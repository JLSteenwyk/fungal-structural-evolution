#!/usr/bin/env python3
"""Independently reconstruct and audit a hydrophobin Pfam candidate inventory."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path


PROFILES = {
    "PF01185.24": ("Hydrophobin", "Family", "chunk-036"),
    "PF06766.17": ("Hydrophobin_2", "Family", "chunk-035"),
    "PF29785.1": ("Hydrophobin_D", "Domain", "chunk-049"),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--domain-db", type=Path, required=True)
    parser.add_argument("--structure-registry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)

    producer = json.loads(args.producer_receipt.read_text())
    if producer["status"] != "complete_provisional_hydrophobin_profile_inventory":
        raise ValueError("Unexpected producer receipt")
    actual_inventory_hash = sha256(args.inventory)
    if producer["artifacts"]["hydrophobin_candidates.tsv"] != actual_inventory_hash:
        raise ValueError("Inventory digest differs from producer receipt")

    catalog = {row["chunk"]: row for row in tsv(args.catalog)}
    if len(catalog) != 64:
        raise ValueError("Full 64-shard catalog is required")
    hits = {}
    source_hashes = {
        str(args.catalog): sha256(args.catalog),
        str(args.domain_db): sha256(args.domain_db),
        str(args.structure_registry): sha256(args.structure_registry),
        str(args.manifest): sha256(args.manifest),
        str(args.inventory): actual_inventory_hash,
        str(args.producer_receipt): sha256(args.producer_receipt),
    }
    for accession, (name, kind, chunk) in PROFILES.items():
        catalog_row = catalog.get(chunk)
        if not catalog_row:
            raise ValueError(f"Catalog lacks {chunk}")
        path = Path(catalog_row["annotation_path"])
        digest = sha256(path)
        if digest != catalog_row["annotation_sha256"]:
            raise ValueError(f"Changed annotation shard {path}")
        source_hashes[str(path)] = digest
        with gzip.open(path, "rt", newline="") as handle:
            for hit in csv.DictReader(handle, delimiter="\t"):
                if hit["pfam_accession"] != accession:
                    continue
                if (hit["pfam_name"], hit["pfam_type"]) != (name, kind):
                    raise ValueError(f"Unexpected Pfam metadata for {accession}")
                if hit["hit_id"] in hits:
                    raise ValueError("Duplicate Pfam hit identity")
                hits[hit["hit_id"]] = hit

    manifest = {row["taxon_id"]: row for row in tsv(args.manifest)}
    if len(manifest) != 526 or {row["study_role"] for row in manifest.values()} != {"ingroup", "outgroup"}:
        raise ValueError("Unexpected frozen analysis manifest")
    records = tsv(args.inventory)
    required = {
        "hit_id", "pfam_accession", "pfam_name", "sequence_id", "taxon_id", "species_name",
        "study_role", "lineage", "protein_id", "query_source", "alignment_start", "alignment_end",
        "envelope_start", "envelope_end", "hmm_coverage", "domain_score", "domain_ga",
        "independent_evalue", "posterior_accuracy", "model_status", "model_id", "model_version",
        "model_length", "model_path",
    }
    if not records or not required.issubset(records[0]):
        raise ValueError("Candidate inventory has required fields missing")
    keys = [(row["hit_id"], row["taxon_id"], row["protein_id"]) for row in records]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate candidate identity")
    if records != sorted(records, key=lambda row: (row["pfam_accession"], row["taxon_id"], row["protein_id"], row["hit_id"])):
        raise ValueError("Candidate order is not canonical")

    domain = sqlite3.connect(f"file:{args.domain_db}?mode=ro", uri=True)
    registry = sqlite3.connect(f"file:{args.structure_registry}?mode=ro", uri=True)
    expected = []
    for hit_id, hit in hits.items():
        protein_rows = domain.execute(
            "SELECT taxon_id, protein_id, query_source FROM proteins WHERE sequence_id=? ORDER BY taxon_id, protein_id",
            (hit["sequence_id"],),
        ).fetchall()
        if not protein_rows:
            raise ValueError(f"Missing protein link: {hit['sequence_id']}")
        model = registry.execute(
            "SELECT model_id, version, length, path FROM models WHERE sequence_sha256=?",
            (hit["sequence_id"].removeprefix("S"),),
        ).fetchone()
        for taxon_id, protein_id, query_source in protein_rows:
            taxon = manifest.get(taxon_id)
            if not taxon:
                raise ValueError(f"Unknown study taxon: {taxon_id}")
            expected.append({
                "hit_id": hit_id, "pfam_accession": hit["pfam_accession"], "pfam_name": hit["pfam_name"],
                "sequence_id": hit["sequence_id"], "taxon_id": taxon_id, "species_name": taxon["species_name"],
                "study_role": taxon["study_role"], "lineage": taxon["lineage"], "protein_id": protein_id,
                "query_source": query_source, "alignment_start": hit["alignment_start"], "alignment_end": hit["alignment_end"],
                "envelope_start": hit["envelope_start"], "envelope_end": hit["envelope_end"],
                "hmm_coverage": hit["hmm_coverage"], "domain_score": hit["domain_score"], "domain_ga": hit["domain_ga"],
                "independent_evalue": hit["independent_evalue"], "posterior_accuracy": hit["posterior_accuracy"],
                "model_status": "source_model_present" if model else "no_source_model",
                "model_id": model[0] if model else "", "model_version": str(model[1]) if model else "",
                "model_length": str(model[2]) if model else "", "model_path": model[3] if model else "",
            })
    expected.sort(key=lambda row: (row["pfam_accession"], row["taxon_id"], row["protein_id"], row["hit_id"]))
    for observed, reconstructed in zip(records, expected):
        for field, value in reconstructed.items():
            if observed[field] != value:
                raise ValueError(f"Candidate differs at {field}: {observed['hit_id']}")
    if len(records) != len(expected):
        raise ValueError("Candidate inventory row count differs from reconstruction")

    profile_counts = dict(sorted(Counter(row["pfam_accession"] for row in records).items()))
    status_counts = dict(sorted(Counter(row["model_status"] for row in records).items()))
    if (len(hits) != producer["accepted_profile_hits"] or len(records) != producer["taxon_protein_profile_rows"]
            or profile_counts != producer["profile_row_counts"] or status_counts != producer["model_status_counts"]):
        raise ValueError("Producer aggregate counts differ")
    result = {
        "status": "passed_independent_hydrophobin_profile_inventory_readback",
        "counts": {"profile_hits": len(hits), "candidate_rows": len(records),
                   "unique_sequences": len({row["sequence_id"] for row in records}),
                   "taxa": len({row["taxon_id"] for row in records})},
        "profile_row_counts": profile_counts,
        "model_status_counts": status_counts,
        "source_hashes": source_hashes,
        "script_sha256": sha256(Path(__file__)),
        "scope": "Every exact hit to the three named profiles was independently rescanned from the catalog shards and every candidate row rebuilt from the frozen manifest, full-domain protein links and source-structure registry. This checks inventory provenance, not homology, structure quality, functional annotation, gene-copy number or evolutionary interpretation.",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
