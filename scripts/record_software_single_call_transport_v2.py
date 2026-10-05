#!/usr/bin/env python3
"""Close an original software call that finished without creating a tool session."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('execution', 'validation', 'payload', 'unit', 'output'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args()
    ep, vp, pp, target = map(Path, (args.execution, args.validation, args.payload, args.output))
    assert not target.exists()
    execution, validation, payload = [json.loads(p.read_text()) for p in (ep, vp, pp)]
    original = payload['original_tool_result']
    assert original['exit_code'] == 0 and 'session_id' not in original and original['chunk_id']
    assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0 and not execution['timed_out']
    assert execution['receipt_sha256'] == sha(vp) and validation['scientific_eligibility'] is False
    invocation = execution['invocation_id']
    assert invocation in original['output'] and args.unit in original['output']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', args.unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [r for r in rows if invocation in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID'))]
    exact = [r for r in rows if r.get('_PID') == str(execution['wrapper']['pid'])
             and r.get('_CMDLINE') == ' '.join(execution['wrapper']['cmdline'])]
    assert len(exact) == 2
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=execution['wrapper'], invocation_id=invocation)
    terminal = {k: v for k, v in execution.items() if k not in ('source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope')}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == invocation and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == invocation and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == 1 and len(ends) in (0, 1)
    assert 'Finished with result: success' in original['output']
    assert 'Main processes terminated with: code=exited/status=0' in original['output']
    assert not any('Failed with result' in r.get('MESSAGE', '') or 'status=1/FAILURE' in r.get('MESSAGE', '') for r in rows)
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True)+'\n')
    pins = {}
    for mapping in (execution['source_hashes'], execution['artifacts'], validation['source_hashes']):
        for path, digest in mapping.items():
            bind(pins, path, digest)
    for path in (ep, vp, pp, journal, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_single_call_software_result_and_whole_wrapper_payloads',
                  checked_utc=datetime.now(timezone.utc).isoformat(), unit=args.unit, invocation_id=invocation,
                  wrapper=execution['wrapper'], original_tool_session_id=None,
                  original_tool_session_created=False, original_tool_chunk_id=original['chunk_id'],
                  original_tool_terminal_exit_code=0, actual_tool_terminal_exit_code=0,
                  validation_sha256=sha(vp), whole_wrapper_initial_and_terminal_payloads_matched=True,
                  manager_start_records=len(starts), manager_completion_records=len(ends),
                  manager_completion_record_available=bool(ends),
                  source_hashes=pins, scientific_eligibility=False,
                  scope='The actual original tool result completed in one call with exit zero and no session. '
                        'Exact whole wrapper payloads, native exit/receipt and invocation manager records '
                        'match. Optional manager resource-completion record availability is explicit; its absence in a fast call is not native failure. Actual original tool exit and whole native wrapper terminal remain required. No session, manager completion record or separate wait is invented; software evidence only.')
    with target.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'wrapper')}, indent=2))


if __name__ == '__main__':
    main()
