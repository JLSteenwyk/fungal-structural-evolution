#!/usr/bin/env python3
"""Run one source-verified taxon through the full-panel RepeatModeler stage."""
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


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(command: list[str], cwd: Path, logfile: Path) -> None:
    with logfile.open("a") as log:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def fasta_census(path: Path) -> tuple[int, int]:
    records = bases = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--build-database", type=Path, required=True)
    parser.add_argument("--repeatmodeler", type=Path, required=True)
    parser.add_argument("--repeatclassifier", type=Path, required=True)
    parser.add_argument("--famdb-dir", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    project_root = config_path.parent.parent
    input_manifest = (project_root / config["input_manifest"]).resolve()
    output_root = (project_root / config["output_root"]).resolve()
    build_database = args.build_database.resolve()
    repeatmodeler = args.repeatmodeler.resolve()
    repeatclassifier = args.repeatclassifier.resolve()
    famdb_dir = args.famdb_dir.resolve()
    if not (famdb_dir / "famdb.py").is_file():
        raise FileNotFoundError(f"FamDB installation unavailable: {famdb_dir}")
    with input_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != config["taxa"] or not 0 <= args.task_index < len(rows):
        raise ValueError("task index or complete-panel manifest is invalid")
    row = rows[args.task_index]
    source = Path(row["genome_fasta_gz"]).resolve()
    if not source.is_file() or digest(source) != row["sha256"]:
        raise ValueError(f"source FASTA verification failed for {row['taxon_id']}")

    work = output_root / row["taxon_id"]
    work.mkdir(parents=True, exist_ok=True)
    receipt = work / "receipt.json"
    if receipt.exists():
        prior = json.loads(receipt.read_text())
        library = Path(prior.get("classified_library", ""))
        if (prior.get("status") == "completed" and prior.get("source_sha256") == row["sha256"]
                and library.is_file() and digest(library) == prior.get("classified_library_sha256")):
            print(f"already completed {row['taxon_id']}")
            return
        raise RuntimeError(f"conflicting completed receipt for {row['taxon_id']}")

    log = work / "run.log"
    fasta = work / "source.fasta"
    if not fasta.exists():
        with gzip.open(source, "rb") as source_handle, fasta.open("wb") as output_handle:
            shutil.copyfileobj(source_handle, output_handle, length=1024 * 1024)
    records, bases = fasta_census(fasta)
    if records != int(row["fasta_records"]) or bases != int(row["fasta_bases"]):
        raise ValueError(f"decompressed FASTA census mismatch for {row['taxon_id']}")

    database = "repeatmodeler_db"
    if not any(work.glob(database + ".*")):
        run([str(build_database), "-name", database, str(fasta)], work, log)
    recover_dirs = sorted(path for path in work.glob("RM_*") if path.is_dir())
    recover_dir = None
    completed_unclassified = None
    if recover_dirs:
        if len(recover_dirs) != 1:
            raise RuntimeError(f"ambiguous recovery directories for {row['taxon_id']}")
        candidate = recover_dirs[0]
        completed_rounds = [
            number for number in range(1, 101)
            if (candidate / f"round-{number}" / "consensi.fa").is_file()
            and (candidate / f"round-{number}" / "consensi.fa").stat().st_size > 0
        ]
        run_log = candidate / "rmod.log"
        if ((candidate / "consensi.fa").is_file() and run_log.is_file()
                and "Program Time:" in run_log.read_text(errors="replace")):
            completed_unclassified = candidate
        elif completed_rounds and max(completed_rounds) > 1:
            recover_dir = candidate
        else:
            archive = work / f"interrupted_{candidate.name}"
            suffix = 1
            while archive.exists():
                archive = work / f"interrupted_{candidate.name}_{suffix}"
                suffix += 1
            os.replace(candidate, archive)
    threads = str(config["threads_per_taxon"])
    seed = str((int(hashlib.sha256(row["taxon_id"].encode()).hexdigest()[:8], 16) % 2147483646) + 1)
    if completed_unclassified:
        command = None
    elif recover_dir:
        command = [str(repeatmodeler), "-recoverDir", str(recover_dir), "-threads", threads, "-srand", seed, "-famdb_dir", str(famdb_dir)]
    else:
        command = [str(repeatmodeler), "-database", database, "-threads", threads, "-srand", seed, "-famdb_dir", str(famdb_dir)]
    if command:
        run(command, work, log)
    candidates = sorted(work.glob("RM_*/consensi.fa.classified"))
    if not candidates:
        raw_candidates = sorted(work.glob("RM_*/consensi.fa"))
        if len(raw_candidates) != 1:
            raise RuntimeError(f"missing completed raw consensus for {row['taxon_id']}")
        classify = [str(repeatclassifier), "-consensi", str(raw_candidates[0]), "-threads", threads, "-famdb_dir", str(famdb_dir)]
        classifier_env = dict(os.environ, FAMDB_DIR=str(famdb_dir))
        with log.open("a") as log_handle:
            log_handle.write("$ " + " ".join(classify) + "\n")
            log_handle.flush()
            subprocess.run(classify, cwd=work, stdout=log_handle, stderr=subprocess.STDOUT, check=True, env=classifier_env)
        candidates = sorted(work.glob("RM_*/consensi.fa.classified"))
    if len(candidates) != 1 or candidates[0].stat().st_size == 0:
        raise RuntimeError(f"missing classified library for {row['taxon_id']}")
    library = candidates[0].resolve()
    payload = {
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
        "repeatmodeler_version": config["repeatmodeler"]["version"],
        "engine": config["repeatmodeler"]["engine"],
        "threads": int(threads),
        "seed": int(seed),
        "classified_library": str(library),
        "classified_library_bytes": library.stat().st_size,
        "classified_library_sha256": digest(library),
        "scope": "Per-taxon de-novo repeat-library discovery. It does not annotate all repeat instances or support genome-architecture inference by itself."
    }
    atomic_json(receipt, payload)
    print(f"completed {row['taxon_id']}")


if __name__ == "__main__":
    main()
