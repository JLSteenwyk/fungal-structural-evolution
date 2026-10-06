#!/usr/bin/env python3
"""Bind the completed coordinate-profile audit to its original systemd terminal."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--initial-tool", type=Path, required=True)
    parser.add_argument("--unit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    execution = json.loads(args.execution.read_text())
    receipt = json.loads(args.receipt.read_text())
    initial = json.loads(args.initial_tool.read_text())
    if execution["status"] != "exited_zero_with_receipt" or execution["exit_code"] != 0 or execution["timed_out"]:
        raise ValueError("Coordinate-profile execution did not close with exit zero")
    if receipt["status"] != "completed_full_atlas_coordinate_and_residue_audit_pending_independent_readback":
        raise ValueError("Unexpected coordinate-profile receipt")
    if execution["receipt_sha256"] != sha(args.receipt):
        raise ValueError("Execution is not bound to coordinate-profile receipt")
    if set(initial) != {"original_tool_session_id", "initial"} or initial["original_tool_session_id"] != initial["initial"]["session_id"]:
        raise ValueError("Unexpected original tool initial record")
    invocation = execution["invocation_id"]
    if invocation not in initial["initial"]["output"] or args.unit not in initial["initial"]["output"]:
        raise ValueError("Original tool initial record does not identify the execution")
    all_rows = [json.loads(line) for line in subprocess.check_output(
        ["journalctl", "--user", "-u", args.unit, "-o", "json", "--no-pager"], text=True).splitlines()]
    rows = [row for row in all_rows if invocation in (row.get("_SYSTEMD_INVOCATION_ID"), row.get("USER_INVOCATION_ID"))]
    exact = [row for row in rows if row.get("_PID") == str(execution["wrapper"]["pid"])
             and row.get("_CMDLINE") == " ".join(execution["wrapper"]["cmdline"])]
    if len(exact) != 2:
        raise ValueError("Expected original wrapper initial and terminal journal rows")
    if json.loads(exact[0]["MESSAGE"]) != {"original_wrapper": execution["wrapper"], "invocation_id": invocation}:
        raise ValueError("Original wrapper initial row differs")
    terminal = {key: value for key, value in execution.items()
                if key not in ("source_hashes", "artifacts", "command", "wrapper", "child", "scope")}
    if json.loads(exact[1]["MESSAGE"]) != terminal:
        raise ValueError("Original wrapper terminal row differs")
    starts = [row for row in rows if row.get("USER_INVOCATION_ID") == invocation and "Started " in row.get("MESSAGE", "")]
    ends = [row for row in rows if row.get("USER_INVOCATION_ID") == invocation and row.get("CPU_USAGE_NSEC")]
    if len(starts) != 1 or len(ends) != 1:
        raise ValueError("Original manager start/terminal rows differ")
    journal = args.execution.with_suffix("") / "original-invocation-journal.jsonl"
    if journal.exists():
        raise FileExistsError(journal)
    with journal.open("x") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    hashes = {}
    # The producer receipt already binds the full source/archive set. Rehashing
    # every entry here would reread the entire atlas merely to close execution
    # provenance, so this closure binds the compact producer receipt instead.
    for mapping in (execution["source_hashes"], execution["artifacts"]):
        for path, digest in mapping.items():
            bind(hashes, path, digest)
    for path in (args.execution, args.receipt, args.initial_tool, journal, Path(__file__)):
        bind(hashes, path)
    verify(hashes)
    result = {
        "status": "verified_original_software_wait_and_whole_wrapper_payloads",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "unit": args.unit,
        "invocation_id": invocation,
        "wrapper": execution["wrapper"],
        "original_tool_session_id": initial["original_tool_session_id"],
        "original_tool_terminal_exit_code": 0,
        "actual_tool_session_id": initial["original_tool_session_id"],
        "actual_tool_terminal_exit_code": 0,
        "validation_sha256": sha(args.receipt),
        "whole_wrapper_initial_and_terminal_payloads_matched": True,
        "manager_start_records": 1,
        "manager_completion_records": 1,
        "exact_wrapper_pid_journal_entries": 2,
        "original_start_records": 1,
        "original_completion_records": 1,
        "source_hashes": hashes,
        "scientific_eligibility": False,
        "scope": "The original external tool initial session record identifies the systemd invocation. "
                 "Its terminal result is established from the original bounded execution receipt and "
                 "matching wrapper/systemd terminal records; no unavailable external wait payload is invented. "
                 "This closes execution provenance only, not independent coordinate readback or biological eligibility.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key not in ("source_hashes", "wrapper")}, indent=2))


if __name__ == "__main__":
    main()
