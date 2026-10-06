#!/usr/bin/env python3
"""Record local repeat-annotation and synteny software availability."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


TOOLS = {
    "RepeatMasker": ["--version"], "RepeatModeler": ["--version"], "BuildDatabase": ["--version"],
    "EDTA": ["--version"], "TEsorter": ["--version"], "Red": ["--version"],
    "MCScanX": ["--version"], "McScanX": ["--version"], "i-adhore": ["--version"],
    "jcvi": ["--version"], "minimap2": ["--version"], "blastp": ["-version"],
    "diamond": ["--version"], "mmseqs": ["version"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    observed = {}
    for tool, version_args in TOOLS.items():
        executable = shutil.which(tool)
        if executable is None:
            observed[tool] = {"status": "not_found"}
            continue
        run = subprocess.run([executable, *version_args], text=True, capture_output=True, timeout=30)
        text = (run.stdout or run.stderr).strip().splitlines()
        observed[tool] = {"status": "available", "path": executable, "sha256": sha256(Path(executable)),
                          "version_command": [executable, *version_args], "version_first_line": text[0] if text else "",
                          "version_exit_code": run.returncode}
    required_missing = [tool for tool in ("RepeatMasker", "RepeatModeler", "EDTA", "MCScanX", "i-adhore", "jcvi")
                        if observed[tool]["status"] == "not_found"]
    result = {
        "schema_version": 1, "status": "completed_genome_architecture_software_availability_audit",
        "checked_utc": datetime.now(timezone.utc).isoformat(), "tools": observed,
        "required_repeat_and_synteny_tools_not_found": required_missing,
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "This audits local executable availability only. Missing repeat/synteny software prevents a full discovery launch until a versioned, no-paid-resource installation or container plan is reviewed; available general aligners do not substitute for repeat annotation or collinearity inference.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"required_repeat_and_synteny_tools_not_found": required_missing}, indent=2))


if __name__ == "__main__":
    main()
