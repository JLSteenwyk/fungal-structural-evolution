#!/usr/bin/env python3
"""Verify the frozen 501-fungus/25-outgroup analysis manifest."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

FUNGAL_LINEAGES = {
    "Aphelidiomycota", "Ascomycota", "Basidiobolomycota", "Basidiomycota",
    "Blastocladiomycota", "Calcarisporiellomycota", "Chytridiomycota",
    "Entomophthoromycota", "Glomeromycota", "Kickxellomycota", "Microsporidia",
    "Monoblepharomycota", "Mortierellomycota", "Mucoromycota",
    "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota",
    "Sanchytriomycota", "Zoopagomycota",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not reader.fieldnames or any(not row.get("taxon_id") for row in rows):
        raise ValueError(f"{path} lacks nonempty taxon_id values")
    return reader.fieldnames, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    candidate_columns, candidate_rows = read_tsv(args.candidate)
    analysis_columns, analysis_rows = read_tsv(args.analysis)
    if candidate_columns != analysis_columns:
        raise ValueError("Candidate and analysis manifest columns differ")
    candidate = {row["taxon_id"]: row for row in candidate_rows}
    analysis = {row["taxon_id"]: row for row in analysis_rows}
    if len(candidate) != len(candidate_rows) or len(analysis) != len(analysis_rows):
        raise ValueError("A manifest has duplicate taxon_id values")
    if not set(analysis) <= set(candidate):
        raise ValueError("Analysis manifest contains a taxon absent from candidate manifest")
    changed = [taxon for taxon in analysis if analysis[taxon] != candidate[taxon]]
    if changed:
        raise ValueError(f"Analysis rows differ from candidate rows: {changed[:5]}")

    roles = Counter(row["study_role"] for row in analysis_rows)
    if len(analysis_rows) != 526 or roles != Counter({"ingroup": 501, "outgroup": 25}):
        raise ValueError(f"Expected 526 = 501 ingroup + 25 outgroup, found {dict(roles)}")
    if set(roles) != {"ingroup", "outgroup"}:
        raise ValueError("Unexpected study role")
    fungal_paths = [row["lineage"].split(";", 1)[0] for row in analysis_rows if row["study_role"] == "ingroup"]
    unknown_lineages = sorted(set(fungal_paths) - FUNGAL_LINEAGES)
    if unknown_lineages:
        raise ValueError(f"Unexpected ingroup lineage prefixes: {unknown_lineages}")
    fungal_species = [row["species_name"].strip() for row in analysis_rows if row["study_role"] == "ingroup"]
    if any(not name for name in fungal_species):
        raise ValueError("An ingroup entry lacks a literal species name")
    duplicated_fungal_species = sorted(name for name, count in Counter(fungal_species).items() if count > 1)
    if duplicated_fungal_species:
        raise ValueError(f"Duplicated literal fungal species names: {duplicated_fungal_species[:5]}")

    exclusions = json.loads(args.exclusions.read_text())
    if not isinstance(exclusions, list):
        raise ValueError("Exclusions must be a list")
    excluded_ids = {row.get("taxon_id") for row in exclusions}
    missing = set(candidate) - set(analysis)
    if missing != excluded_ids or len(missing) != 1:
        raise ValueError("Candidate-minus-analysis rows do not match the documented exclusion")
    excluded_id = next(iter(missing))
    exclusion = next(row for row in exclusions if row["taxon_id"] == excluded_id)
    for field in ("species_name", "assembly_accession"):
        if exclusion.get(field) != candidate[excluded_id].get(field):
            raise ValueError(f"Exclusion {field} does not match candidate manifest")
    if not exclusion.get("reason") or not exclusion.get("evidence"):
        raise ValueError("Exclusion lacks reason or evidence")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "passed_frozen_analysis_manifest_composition_and_exclusion_audit",
        "analysis_rows": len(analysis_rows),
        "analysis_role_counts": dict(sorted(roles.items())),
        "candidate_rows": len(candidate_rows),
        "candidate_role_counts": dict(sorted(Counter(row["study_role"] for row in candidate_rows).items())),
        "analysis_identity_matches_candidate_rows": True,
        "ingroup_literal_species_names": len(fungal_species),
        "ingroup_literal_species_names_unique": True,
        "excluded_candidate": {
            "taxon_id": excluded_id,
            "species_name": candidate[excluded_id]["species_name"],
            "reason": exclusion["reason"],
            "evidence": exclusion["evidence"],
        },
        "ingroup_top_level_lineage_counts": dict(sorted(Counter(fungal_paths).items())),
        "allowed_fungal_top_level_lineages": sorted(FUNGAL_LINEAGES),
        "source_hashes": {
            str(args.candidate): sha256(args.candidate),
            str(args.analysis): sha256(args.analysis),
            str(args.exclusions): sha256(args.exclusions),
        },
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "The analysis universe contains exactly 501 fungal ingroup entries and 25 outgroups. This validates manifest composition and the one documented candidate exclusion; it does not validate taxonomy, assembly quality, annotation quality, ecological labels, species delimitation, or suitability for each downstream family analysis.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
