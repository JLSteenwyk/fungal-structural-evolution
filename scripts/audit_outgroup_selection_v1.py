#!/usr/bin/env python3
"""Audit frozen outgroup membership against its configuration and manifest.

This checks the reproducible *selection record*.  It deliberately does not
infer a phylogeny or use the configured role labels as phylogenetic evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


EXPECTED_GROUP_COUNTS = {
    "close_holomycotan_anchor": 1,
    "unicellular_holozoan_comparisons": 9,
    "animal_lineage_coverage": 9,
    "broader_rooting_sensitivity": 6,
}
GROUP_BY_LINEAGE = {
    "Nucleariida": "close_holomycotan_anchor",
    "Choanoflagellatea": "unicellular_holozoan_comparisons",
    "Filasterea": "unicellular_holozoan_comparisons",
    "Ichthyosporea": "unicellular_holozoan_comparisons",
    "Corallochytrea": "unicellular_holozoan_comparisons",
    "Metazoa": "animal_lineage_coverage",
    "Apusomonadida": "broader_rooting_sensitivity",
    "Amoebozoa": "broader_rooting_sensitivity",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with args.config.open() as handle:
        config = json.load(handle)

    outgroups = [row for row in rows if row["study_role"] == "outgroup"]
    if len(outgroups) != 25:
        raise ValueError(f"expected exactly 25 manifest outgroups, found {len(outgroups)}")
    names = [row["species_name"] for row in outgroups]
    if len(set(names)) != len(names):
        raise ValueError("manifest outgroup species names are not unique")

    configured = list(config["ncbi"]) + list(config["external"])
    config_names = [entry[0] if isinstance(entry[0], str) else entry[1] for entry in configured]
    if len(config_names) != 25 or len(set(config_names)) != 25:
        raise ValueError("configuration must contain 25 unique selected species names")
    if set(names) != set(config_names):
        missing = sorted(set(config_names) - set(names))
        unexpected = sorted(set(names) - set(config_names))
        raise ValueError(f"configuration/manifest membership mismatch: missing={missing}; unexpected={unexpected}")

    configured_lineage = {}
    for entry in config["ncbi"]:
        name, lineage, _reason = entry
        configured_lineage[name] = lineage
    for entry in config["external"]:
        _figshare_id, name, lineage = entry
        configured_lineage[name] = lineage

    membership = []
    for row in sorted(outgroups, key=lambda item: item["species_name"]):
        top_level = row["lineage"].split(";", 1)[0]
        if top_level not in GROUP_BY_LINEAGE:
            raise ValueError(f"unmapped top-level lineage for {row['species_name']}: {top_level}")
        if configured_lineage[row["species_name"]].split(";", 1)[0] != top_level:
            raise ValueError(f"lineage mismatch for {row['species_name']}")
        membership.append({
            "taxon_id": row["taxon_id"],
            "species_name": row["species_name"],
            "top_level_lineage": top_level,
            "selection_group": GROUP_BY_LINEAGE[top_level],
            "assembly_accession": row["assembly_accession"],
            "inclusion_reason": row["inclusion_reason"],
        })

    group_counts = dict(sorted(Counter(item["selection_group"] for item in membership).items()))
    if group_counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(f"unexpected group counts: {group_counts}")
    evidence = config.get("evidence", [])
    if len(evidence) < 2 or any(not url.startswith("https://doi.org/") for url in evidence):
        raise ValueError("configuration must retain at least two DOI evidence URLs")

    result = {
        "schema_version": 1,
        "audit": "frozen_outgroup_selection_membership_and_provenance",
        "manifest_outgroup_count": len(outgroups),
        "configured_outgroup_count": len(config_names),
        "selection_group_counts": group_counts,
        "top_level_lineage_counts": dict(sorted(Counter(item["top_level_lineage"] for item in membership).items())),
        "membership": membership,
        "evidence_urls": evidence,
        "source_hashes": {
            str(args.manifest): sha256(args.manifest),
            str(args.config): sha256(args.config),
        },
        "script_sha256": sha256(Path(__file__)),
        "interpretation": (
            "This receipt validates the exact frozen outgroup membership, configured lineage labels, "
            "role balance, and cited-source URLs. It does not independently establish taxonomy, "
            "assembly/annotation quality, phylogenetic topology, root placement, divergence times, "
            "or suitability of every outgroup for every protein family."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("manifest_outgroup_count", "selection_group_counts", "top_level_lineage_counts")}, indent=2))


if __name__ == "__main__":
    main()
