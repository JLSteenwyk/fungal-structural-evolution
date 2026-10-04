#!/usr/bin/env python3
"""Bind actual software refits to the original tool wait and manager journal."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--actual-tool-session', type=int, required=True)
    p.add_argument('--actual-tool-exit', type=int, required=True)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--validation', type=Path, required=True)
    p.add_argument('--unit', required=True);p.add_argument('--output', type=Path, required=True)
    a = p.parse_args();assert a.actual_tool_exit == 0
    e = json.loads(a.execution.read_text());v = json.loads(a.validation.read_text())
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and not e['timed_out']
    assert e['receipt_sha256'] == sha(a.validation)
    assert v['actual_software_refits'] == v['serialized_independent_refit_replays'] == 128
    assert v['production_biological_fits_computed'] == 0 and v['scientific_eligibility'] is False
    assert sum(v['refit_status_counts'].values()) == 128
    rows = [json.loads(l) for l in subprocess.check_output(
        ['journalctl', '--user', '-u', a.unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    inv = e['invocation_id']
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID')]]
    exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid'])
        and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(r['_SYSTEMD_INVOCATION_ID'] == inv for r in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: value for k, value in e.items()
        if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = a.execution.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as f:
        for r in rows:f.write(json.dumps(r, sort_keys=True)+'\n')
    bindings = {}
    for mapping in [e['source_hashes'], e['artifacts'], v['source_hashes']]:
        for name, d in mapping.items():bind(bindings, name, d)
    for path in [a.execution, a.validation, journal, Path(__file__)]:bind(bindings, path)
    verify(bindings)
    result = dict(status='verified_actual_weighted_gaussian_refits_original_software_wait_zero',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_tool_session_id=a.actual_tool_session,
        actual_tool_terminal_exit_code=a.actual_tool_exit, validation_sha256=sha(a.validation),
        unit=a.unit, invocation_id=inv, wrapper=e['wrapper'], exact_wrapper_pid_journal_entries=2,
        original_start_records=1, original_completion_records=1, entire_terminal_payload_matched=True,
        source_hashes=bindings, actual_software_refits=128, production_biological_fits_computed=0,
        scientific_eligibility=False,
        scope='Actual original bounded software wait, exact wrapper PID/create/command/invocation, '
            'entire original terminal receipt/hash payload and manager start/completion. '
            'Native mathematics on synthetic fixtures is not full production fitting or calibration.')
    atomic(a.output, result)
    print(json.dumps({k: value for k, value in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__':main()
