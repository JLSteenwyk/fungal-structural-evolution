#!/usr/bin/env python3
"""Preserve the original rejected V1 checker and qualified V2 dispatch controls."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', choices=['v1', 'v2'], required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    ep = Path('metadata/retained_factor_v9_dispatch_software_execution_20261004_'+a.version+'.json')
    vp = Path('metadata/retained_factor_v9_dispatch_software_validation_20261004_'+a.version+'.json')
    e = json.loads(ep.read_text()); successful = a.version == 'v2'
    assert not e['timed_out'] and e['exit_code'] == (0 if successful else 1)
    assert e['status'] == ('exited_zero_with_receipt' if successful else 'failed_stage_retained')
    pins = {}
    for mapping in [e['source_hashes'], e['artifacts']]:
        for name, digest in mapping.items(): bind(pins, name, digest)
    if successful:
        v = json.loads(vp.read_text()); assert sha(vp) == e['receipt_sha256']
        assert v['status'] == 'passed_actual_native_v9_dispatch_metadata_controls_v2'
        assert v['actual_native_svd_calls'] == 3 and v['original_probe_and_guard_runs'] == 0
        assert v['artificial_probe_entry'] and not v['scientific_eligibility']
        assert [row['case'] for row in v['cases']] == ['c_order', 'f_order', 'source_alias']
        assert all(row['actual_native_svd_calls'] == 1 and row['dispatch_line_events'] == 2
            and row['inputs_preserved'] for row in v['cases'])
        for name, digest in v['source_hashes'].items(): bind(pins, name, digest)
        bind(pins, vp)
    else:
        assert not vp.exists() and e['receipt_sha256'] is None
        stderr = ep.with_suffix('')/'stderr.log'
        assert 'assert len(dispatches) == 1' in stderr.read_text()
        root = Path('data/software_audits/retained-factor-v9-dispatch-20261004-v1')
        assert {q.name for q in root.iterdir()} == {'c_order.jsonl'}
        trace = root/'c_order.jsonl'; rows = [json.loads(line) for line in trace.read_text().splitlines()]
        dispatch = [row['native_svd_dispatch'] for row in rows if 'native_svd_dispatch' in row]
        assert len(dispatch) == 2 and dispatch[0] == dispatch[1]
        assert all(not row['changed'] for row in rows)
        bind(pins, trace)
    unit = 'fungal-retained-v9-dispatch-software-20261004-'+a.version+'.service'
    inv = e['invocation_id']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [row for row in rows if inv in [row.get('_SYSTEMD_INVOCATION_ID'), row.get('USER_INVOCATION_ID')]]
    exact = [row for row in rows if row.get('_PID') == str(e['wrapper']['pid'])
        and row.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2 and all(row['_SYSTEMD_INVOCATION_ID'] == inv for row in exact)
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: val for k, val in e.items() if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and 'Started ' in row.get('MESSAGE', '')]
    ends = [row for row in rows if row.get('USER_INVOCATION_ID') == inv and row.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in rows: handle.write(json.dumps(row, sort_keys=True)+'\n')
    for path in [Path(__file__), ep, journal]: bind(pins, path)
    verify(pins)
    payloads = {
        'v1': dict(chunk_id='e518c9', wall_time_seconds=0.000013862, exit_code=1, original_token_count=44,
            output='Finished with result: exit-code\nMain processes terminated with: code=exited/status=1\n'
                   'Service runtime: 1.233s\nCPU time consumed: 1.194s\nMemory peak: 256.0K\nMemory swap peak: 0B\n'),
        'v2': dict(chunk_id='ca002a', wall_time_seconds=0.000018628, exit_code=0, original_token_count=44,
            output='Finished with result: success\nMain processes terminated with: code=exited/status=0\n'
                   'Service runtime: 1.187s\nCPU time consumed: 1.147s\nMemory peak: 256.0K\nMemory swap peak: 0B\n')}
    result = dict(status='verified_original_v9_dispatch_software_' + ('wait_zero' if successful else 'checker_rejection_retained'),
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_tool_session_id=39674 if successful else 55345,
        actual_tool_terminal_exit_code=0 if successful else 1, unit=unit, invocation_id=inv,
        wrapper=e['wrapper'], original_tool_terminal_payload=payloads[a.version],
        entire_terminal_payload_matched=True, original_start_records=1, original_completion_records=1,
        source_hashes=pins, qualified_actual_native_svd_calls=3 if successful else 0,
        recorded_repeated_dispatch_line_events=True, original_factor_runs=0,
        original_jobs_restarted=False, scientific_eligibility=False, biological_fits=0, gpu=False,
        scope='Exact original capped dispatch checker wait and whole wrapper journal payloads plus manager '
              'start/end preserved. V1 rejected repeated CPython pre-dispatch line events; actual SVD '
              'observer stays unchanged. V2 verifies all three actual small SVD calls, including source '
              'alias, with duplicated metadata explicitly retained rather than counted as more calls. '
              'Private substituted probe entry, not original guards, factor corruption, repair or biology.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: val for k, val in result.items() if k not in ['source_hashes', 'wrapper']}, indent=2))


if __name__ == '__main__': main()
