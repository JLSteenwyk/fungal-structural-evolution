#!/usr/bin/env python3
"""Rehash the complete numerical archive and close original reader/closer waits."""
import ast
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import bind, verify


def terminal(label, expected_summary, native=None):
    launch_path = Path(f'metadata/full_weighted_parallel_{label}_launch_20261004_v1.json')
    tool_path = Path(f'metadata/full_weighted_parallel_{label}_original_tool_payloads_20261004_v1.json')
    launch, tool = [json.loads(p.read_text()) for p in [launch_path, tool_path]]
    launch['launch'] = str(launch_path)
    assert fingerprint(launch) is None
    assert tool['original_tool_session_id'] == tool['initial']['session_id'] == launch['original_tool_session_id']
    assert tool['terminal']['exit_code'] == 0
    assert 'Finished with result: success' in tool['terminal']['output']
    assert launch['invocation_id'] in tool['initial']['output']
    assert sha(launch['plan']) == launch['plan_sha256']
    observed = journal_terminal(launch)
    raw = subprocess.check_output(['journalctl', '--user', '-u', launch['unit'],
                                   '--all', '-o', 'json', '--no-pager'], text=True)
    inv = launch['invocation_id']
    rows = [json.loads(line) for line in raw.splitlines()]
    rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID')]]
    wrapper = [r for r in rows if r.get('_PID') == str(launch['pid'])
               and r.get('_CMDLINE') == ' '.join(launch['cmdline'])]
    headers = [r['MESSAGE'] for r in wrapper
               if isinstance(r.get('MESSAGE'), str) and r['MESSAGE'].startswith('running_original_journal_verified_command ')]
    assert len(headers) == 1
    wait = json.loads(Path(launch['plan']).read_text())
    command = ast.literal_eval(headers[0].split('running_original_journal_verified_command ', 1)[1])
    assert command == wait['command']
    expected_native_command = command[command.index('--') + 1:]
    native_rows = [r for r in rows if r.get('_CMDLINE') == ' '.join(expected_native_command)]
    assert native_rows and len({r['_PID'] for r in native_rows}) == 1
    if native is not None:
        assert native['cmdline'] == expected_native_command
        assert {r['_PID'] for r in native_rows} == {str(native['pid'])}
    messages = [r['MESSAGE'] for r in native_rows
                if not r['MESSAGE'].startswith('parallel_weighted_covariance_cohorts ')]
    assert messages
    assert all(r.get('_LINE_BREAK') == 'line-max' for r in native_rows[-len(messages):-1])
    assert json.loads(''.join(messages)) == expected_summary
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv
              and (r.get('MESSAGE') or '').startswith('Started ')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    journal = Path(f'metadata/full_weighted_parallel_{label}_original_whole_journal_20261004_v1.jsonl')
    with journal.open('x') as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + '\n')
    return dict(stage=label, original_tool_session_id=launch['original_tool_session_id'],
                original_tool_terminal_exit_code=0, invocation_id=inv,
                wrapper_pid=launch['pid'], native_pid=int(native_rows[0]['_PID']),
                native_summary_chunks=len(messages), whole_native_summary_matched=True,
                manager_start_records=1, manager_completion_records=1,
                terminal_handle=observed, journal=str(journal), journal_sha256=sha(journal)), [launch_path, tool_path, Path(launch['plan']), journal]


def main():
    output = Path('metadata/full_weighted_parallel_qualification_full_closure_verified_20261004_v1.json')
    assert not output.exists()
    completed_path = Path('metadata/full_weighted_parallel_qualification_completed_20261004_v1.json')
    completed = json.loads(completed_path.read_text())
    archive_path = Path(completed['full_hash_archive'])
    assert sha(archive_path) == completed['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    assert len(archive['source_hashes']) == completed['bound_source_hashes'] == 30561
    assert len(archive['services']) == completed['exact_process_journals_checked'] == 2
    verify(archive['source_hashes'])
    assert completed['cohorts'] == 4340 and completed['designs'] == 130200
    assert completed['numerical_audit_rows'] == 5208000 and completed['setting_audit_links'] == 24883200
    status = 'numerically_qualified_exact_four_control_covariance_basis'
    assert completed['audit_status_counts'] == {status: 5208000}
    assert completed['link_status_counts'] == {status: 24883200}
    assert completed['policies'] == ['uniform', 'background_node', 'background_pair', 'family_component']
    assert completed['loading_modes'] == ['signed', 'unsigned'] and len(set(completed['trees'])) == 5
    assert completed['policy_audit_status_counts'] == {p + ':' + status: 1302000 for p in completed['policies']}
    assert completed['cached_numeric_inputs_preserved_for_all_cohorts'] is True
    assert completed['scientific_eligibility'] is archive['scientific_eligibility'] is False
    root = Path(completed['producer_receipt']).parent
    producer, reader = [json.loads(Path(completed[k]).read_text()) for k in ['producer_receipt', 'independent_readback']]
    for key in ['producer_receipt', 'independent_readback']:
        assert sha(completed[key]) == completed[key + '_sha256']
    assert reader['producer_receipt_sha256'] == sha(completed['producer_receipt'])
    assert producer['working_model_fits_computed'] == reader['working_model_fits_computed'] == 0
    assert reader['conservative_diagnostic_disagreements'] == 0
    for field, value in archive['summary'].items():
        assert producer[field] == reader[field] == completed[field] == value
    assert len(reader['reader_checkpoint_artifacts']) == 4340
    worker_ids = {}
    for label in ['producer', 'reader']:
        checkpoints = sorted((root / 'checkpoints').glob('*.' + label + '.json'))
        assert len(checkpoints) == 4340
        identities = set()
        for index, path in enumerate(checkpoints):
            cp = json.loads(path.read_text())
            assert cp['cohort_index'] == index and cp['cached_numeric_inputs_preserved']
            assert cp['cached_input_array_bindings_checked'] > 0
            assert cp['worker_address_space_limit_bytes'] == 12 * 2**30
            identities.add((cp['worker']['pid'], cp['worker']['created']))
        assert len(identities) == 16
        worker_ids[label] = sorted(identities)
    producer_proof = Path('metadata/full_weighted_parallel_producer_terminal_verified_20261004_v2.json')
    prior = json.loads(producer_proof.read_text()); verify(prior['source_hashes'])
    assert prior['original_tool_terminal_exit_code'] == 0
    assert prior['whole_original_command_header_and_native_terminal_summary_matched']
    snapshot = Path('metadata/full_weighted_parallel_reader_checkpoint_20261004_goal_1306.json')
    native = json.loads(snapshot.read_text())['native_parent']
    reader_summary = {k: v for k, v in reader.items() if k not in ['source_hashes', 'artifacts', 'scope']}
    closures = []; paths = []
    for label, summary, identity in [('reader', reader_summary, native), ('closure', {k: v for k, v in completed.items() if k != 'scope'}, None)]:
        proof, bound = terminal(label, summary, identity); closures.append(proof); paths.extend(bound)
    pins = {}
    for path in [Path(__file__), completed_path, archive_path, producer_proof, snapshot, *paths]:
        bind(pins, path)
    verify(pins)
    result = dict(status='fresh_complete_full_parallel_numerical_archive_and_original_waits_verified',
                  checked_utc=datetime.now(timezone.utc).isoformat(), archive_bindings_freshly_rehashed=30561,
                  producer_bindings_freshly_rehashed=len(prior['source_hashes']), cohorts=4340,
                  numerical_audit_rows=5208000, setting_audit_links=24883200,
                  reader_checkpoint_count=4340, original_worker_identities=worker_ids,
                  original_reader_and_closer_transport=closures,
                  full_numerical_qualification_closed=True, source_hashes=pins,
                  biological_fits=0, scientific_eligibility=False, all_eight_aims_incomplete=True,
                  original_jobs_restarted=False, gpu=False,
                  scope='Complete archive/source/output rehash and checkpoint census, original actual waits, '
                        'wrapper headers, complete distinct-native terminal summaries (including nine reader '
                        'journal chunks), and invocation-bound manager records. Original independent reader '
                        'performed full arithmetic; this check does not recompute it, repair historical SVD '
                        'mutation, optimize biological models, accept nonuniform inference, or qualify posterior uncertainty.')
    with output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_reader_and_closer_transport', 'original_worker_identities']}, indent=2))


if __name__ == '__main__':
    main()
