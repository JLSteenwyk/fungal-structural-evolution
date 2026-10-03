#!/usr/bin/env python3
"""Observe closed startup/diagnostic stages and original software invocations.

This checks compact archive bindings and immutable plan pins. It does not
repeat the completed full native/diagnostic readers or hash every large native
artifact again, and it cannot qualify an ancestral posterior.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_reference_preflight import verify


STAGES = [
    ("reference_preflight", "baliphy_reference_preflight", "20261003"),
    ("reference_footer", "baliphy_reference_startup_footer", "20261003"),
    ("categorical", "independent_baliphy_category_v2", "20261002"),
]
SOFTWARE = [
    ("fungal-sampler-resource-software-validation-20261003-v2.service",
     "metadata/baliphy_sampler_resource_observation_software_validation_20261003_v2.json",
     "scripts/check_baliphy_sampler_resource_observation.py"),
    ("fungal-sampler-resource-serialization-validation-20261003-v1.service",
     "metadata/baliphy_sampler_resource_observation_combined_validation_20261003.json",
     "scripts/check_baliphy_sampler_resource_serialization.py"),
    ("fungal-completed-category-review-tables-20261003-v2.service",
     "metadata/ancestral_category_completed_review_tables_20261003.json",
     "scripts/export_completed_ancestral_category_reviews.py"),
]


def software_terminal(unit, receipt_path, script):
    """Deduplicate journal identities; require actual invocation-linked resources."""
    raw = subprocess.check_output(
        ["journalctl", "--user", "-u", unit, "-o", "json", "--no-pager"], text=True)
    rows = [json.loads(line) for line in raw.splitlines()]
    identities = sorted({(r["_PID"], r["_CMDLINE"], r["_SYSTEMD_INVOCATION_ID"])
                         for r in rows if script in r.get("_CMDLINE", "")
                         and r.get("_SYSTEMD_INVOCATION_ID")})
    assert len(identities) == 1, (unit, identities)
    pid, command, invocation = identities[0]
    # A collected unit is not terminal proof. Original process absence plus
    # the matching manager completion/resource record and finished receipt are.
    assert not Path("/proc", pid).exists(), (unit, "Original PID needs live/reuse review")
    resources = [r for r in rows if r.get("USER_INVOCATION_ID") == invocation
                 and r.get("CPU_USAGE_NSEC")]
    failures = [r for r in rows if r.get("USER_INVOCATION_ID") == invocation
                and ("Main process exited" in r.get("MESSAGE", "")
                     or "Failed with result" in r.get("MESSAGE", ""))]
    assert resources and not failures, unit
    receipt = json.loads(Path(receipt_path).read_text())
    assert receipt["status"].startswith(("passed_", "exported_full_closed_"))
    assert receipt["scientific_eligibility"] is False
    verify({"pins": receipt["source_hashes"]})
    for path, digest in receipt.get("artifacts", {}).items():
        assert sha(path) == digest, path
    return dict(unit=unit, pid=int(pid), command=command, invocation=invocation,
                receipt=receipt_path, receipt_sha256=sha(receipt_path),
                checked_source_bindings=len(receipt["source_hashes"]),
                checked_export_bindings=len(receipt.get("artifacts", {})),
                journal_sha256=__import__("hashlib").sha256(raw.encode()).hexdigest(),
                actual_completion_resources=[{k: r.get(k) for k in
                    ["__REALTIME_TIMESTAMP", "CPU_USAGE_NSEC", "MEMORY_PEAK", "MEMORY_SWAP_PEAK"]}
                    for r in resources], status="verified_original_software_completion")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    stages = {}
    for name, stem, date in STAGES:
        plan_path = Path(f"metadata/{stem}_plan_{date}.json")
        inventory_path = Path(f"metadata/{stem}_launches_{date}.json")
        completed_path = Path(f"metadata/{stem}_completed_{date}.json")
        plan = json.loads(plan_path.read_text())
        verify(plan)
        inventory = json.loads(inventory_path.read_text())
        assert inventory["source_plan_sha256"] == sha(plan_path)
        completed = json.loads(completed_path.read_text())
        archive_path = Path(completed["full_hash_archive"])
        assert sha(archive_path) == completed["full_hash_archive_sha256"]
        archive = json.loads(archive_path.read_text())
        assert len(archive["source_hashes"]) == completed["bound_source_hashes"]
        assert len(archive["services"]) == completed["exact_process_journals_checked"] == 2
        assert completed["scientific_eligibility"] is archive["scientific_eligibility"] is False
        for key, value in archive["summary"].items():
            assert completed[key] == value, (name, key)
        for key in ["producer_receipt", "independent_readback"]:
            path = completed[key]
            assert sha(path) == completed[key + "_sha256"]
            assert archive["source_hashes"][path] == completed[key + "_sha256"]
        handles = []
        for launch in inventory["launches"]:
            record = json.loads(Path(launch).read_text())
            record["launch"] = launch
            assert sha(record["plan"]) == record["plan_sha256"]
            assert fingerprint(record) is None, (name, "Original handle remains live")
            handles.append(journal_terminal(record))
        assert len(handles) == 3
        assert {service["launch"] for service in archive["services"]} == set(inventory["launches"][:2])
        stages[name] = dict(plan=str(plan_path), plan_sha256=sha(plan_path),
            frozen_plan_pins_checked=len(plan["pins"]), completion=str(completed_path),
            completion_sha256=sha(completed_path), archive=str(archive_path),
            archive_sha256=sha(archive_path), archive_binding_count=len(archive["source_hashes"]),
            summary=archive["summary"], original_terminal_handles=handles)
    assert stages["reference_footer"]["summary"]["validated_startups"] == 1620
    assert stages["reference_footer"]["summary"]["unsuccessful_startups"] == 0
    assert stages["reference_footer"]["summary"]["native_timing_footers_checked"] == 72
    assert stages["categorical"]["summary"]["full_quartets"] == 405
    assert stages["categorical"]["summary"]["binary_zero_pair_indicator_rows"] == 15
    result = dict(status="verified_closed_ancestral_qualification_checkpoint",
        checked_utc=datetime.now(timezone.utc).isoformat(), observer_sha256=sha(__file__),
        stages=stages, software_invocations=[software_terminal(*spec) for spec in SOFTWARE],
        posterior_qualified=False, all_eight_aims_incomplete=True, new_native_jobs=0, gpu=False,
        scope="All frozen plan pins, compact completion archive SHA and summaries, producer/reader receipt bindings and all nine original terminal handles checked. Original full archive closures remain the evidence for large native artifacts; this observer does not rehash them or rerun full readers. Three original software/export invocations have deduplicated actual journal identities, source/export hashes and terminal resources. No restart, allocation repair, posterior or biological acceptance.")
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ["stages", "software_invocations"]}))


if __name__ == "__main__":
    main()
