#!/usr/bin/env python3
"""Close original native extraction while retaining its unavailable API terminal."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


NAMES = dict(launch='metadata/completed_domain_coordinates_launch_20261005_v1.json',
    initial='metadata/completed_domain_coordinates_original_initial_tool_20261005_v1.json',
    unavailable='metadata/completed_domain_coordinates_original_api_unavailable_20261006_v1.json',
    execution='metadata/completed_domain_coordinates_execution_20261005_v1.json',
    producer='metadata/completed_domain_coordinates_20261005_v1.json',
    resources='metadata/completed_domain_coordinates_resources_20261005_v1.json',
    plan='metadata/completed_domain_coordinates_plan_20261005_v1.json')


def inputs():
    docs = {name: json.loads(Path(path).read_text()) for name, path in NAMES.items()}
    ep = Path(NAMES['execution'])
    docs['configuration'] = json.loads((ep.with_suffix('') / 'configuration.json').read_text())
    docs['process'] = json.loads((ep.with_suffix('') / 'process.json').read_text())
    docs['raw'] = json.loads((Path(docs['plan']['output']) / 'receipt.json').read_text())
    unit, inv = docs['launch']['unit'], docs['launch']['invocation_id']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True).splitlines()]
    docs['journal'] = [r for r in rows if inv in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID'))]
    return docs


def validate(docs):
    e, v, launch, initial, unavailable, cfg, child, raw, resources, plan = [docs[k] for k in
        ['execution', 'producer', 'launch', 'initial', 'unavailable', 'configuration', 'process', 'raw', 'resources', 'plan']]
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and e['timed_out'] is False
    assert v['status'] == 'completed_full_refreshed_domain_coordinates_pending_independent_atom_readback'
    assert v['scientific_eligibility'] is False
    assert e['receipt'] == NAMES['producer'] and e['receipt_sha256'] == sha(NAMES['producer'])
    assert set(initial) == {'original_tool_session_id', 'initial'}
    assert initial['initial']['session_id'] == initial['original_tool_session_id'] == launch['original_tool_session_id'] == 64749
    assert unavailable['original_tool_session_id'] == 64749 and unavailable['tool_call_failed'] is True
    assert 'Unknown process id 64749' in unavailable['error'] and 'exit_code' not in unavailable
    assert launch['invocation_id'] == cfg['invocation_id'] == e['invocation_id']
    assert launch['pid'] == cfg['wrapper']['pid'] == e['wrapper']['pid']
    assert launch['created'] == cfg['wrapper']['created'] == e['wrapper']['created']
    assert launch['cmdline'] == cfg['wrapper']['cmdline'] == e['wrapper']['cmdline']
    assert launch['unit'] in initial['initial']['output'] and e['invocation_id'] in initial['initial']['output']
    assert child == e['child'] and child['command'] == cfg['command'] == e['command']
    for key in ('working_directory', 'timeout_seconds', 'started_utc', 'source_hashes', 'actual_cgroup_limits'):
        assert cfg[key] == e[key]
    caps = {'cpu.max': '400000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
    assert e['actual_cgroup_limits'] == launch['actual_cgroup_limits'] == caps
    assert resources['cpus'] == 4 and resources['memory_gib'] == 32 and resources['swap_gib'] == 0
    assert e['timeout_seconds'] == resources['wall_seconds_per_stage'] == 345600
    command = e['command']
    assert command[:5] == ['/usr/bin/prlimit', '--as=' + str(resources['address_space_gib'] * 2**30),
        '--cpu=' + str(resources['cpu_seconds_per_stage']), '--fsize=' + str(resources['per_file_limit_mib'] * 2**20), '--']
    assert command[6:] == ['scripts/run_completed_domain_coordinates_v1.py', '--plan', NAMES['plan'], '--receipt', NAMES['producer']]
    rows, inv = docs['journal'], e['invocation_id']
    assert rows and all(inv in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID')) for r in rows)
    exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid']) and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
    assert len(exact) == 2
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
    terminal = {k: v for k, v in e.items() if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
    assert json.loads(exact[1]['MESSAGE']) == terminal
    starts = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and 'Started ' in r.get('MESSAGE', '')]
    ends = [r for r in rows if r.get('USER_INVOCATION_ID') == inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts) == len(ends) == 1
    assert int(ends[0]['CPU_USAGE_NSEC']) > 0 and int(ends[0]['__REALTIME_TIMESTAMP']) > int(starts[0]['__REALTIME_TIMESTAMP'])
    assert raw['status'] == 'complete_domain_coordinate_dispositions_pending_independent_readback'
    assert raw['plan_sha256'] == launch['plan_sha256'] == sha(NAMES['plan'])
    assert sha(Path(plan['output']) / 'receipt.json') == v['raw_receipt_sha256']
    assert raw['counts'] == v['counts'] and raw['shards'] == v['shards'] == len(raw['proofs']) == 977
    assert raw['counts']['models'] == plan['expected_models'] == 976357
    assert raw['counts']['intervals'] == plan['expected_intervals'] == 2454565
    assert raw['counts']['exported'] + raw['counts']['rejected'] == 2454565
    aggregate, jobs = Counter(), set()
    for proof in raw['proofs']:
        assert proof['job'] not in jobs
        jobs.add(proof['job'])
        aggregate.update(proof['receipt']['counts'])
        assert all(isinstance(h, str) and len(h) == 64 for h in proof['receipt']['artifacts'].values())
    assert all(aggregate[k] == v for k, v in raw['counts'].items())
    return dict(original_tool_session_id=64749, original_tool_terminal_exit_code=None,
        original_tool_terminal_available=False, original_native_exit_code=0,
        original_native_exit_proof_available=True, whole_wrapper_initial_and_terminal_payloads_matched=True,
        manager_start_records=1, manager_completion_records=1, exact_wrapper_pid_journal_entries=2,
        counts=raw['counts'], shards=977, native_work_repeated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qualification', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    qualified = json.loads(args.qualification.read_text())
    assert qualified['status'] == 'passed_domain_native_closure_actual_records_and_identity_rejections'
    assert qualified['closure_script_sha256'] == sha(Path(__file__))
    verify(qualified['source_hashes'])
    docs = inputs()
    checked = validate(docs)
    pins = dict(qualified['source_hashes'])
    e, v = docs['execution'], docs['producer']
    for mapping in (e['source_hashes'], e['artifacts']):
        for p, h in mapping.items():
            assert sha(p) == h
            bind(pins, p, h)
    root = Path('results/domains/full-domain-native-closure-20261006-v1')
    root.mkdir(parents=True, exist_ok=False)
    journal = root / 'original-invocation-journal.jsonl'
    with journal.open('x') as handle:
        for row in docs['journal']:
            handle.write(json.dumps(row, sort_keys=True) + '\n')
    archive = root / 'producer-source-bindings.json'
    with archive.open('x') as handle:
        json.dump(dict(status='original_producer_verified_source_bindings_not_fresh_coordinate_readback',
            source_hashes=v['source_hashes']), handle, indent=2)
        handle.write('\n')
    # The original producer already rehashed all declared bulk files before its
    # terminal receipt. Avoid another bulk pass here; independent atom readback
    # will freshly hash every shard and used CIF. This distinction is explicit.
    for p in [*map(Path, NAMES.values()), args.qualification, journal, archive, Path(__file__),
              Path(docs['plan']['output']) / 'receipt.json']:
        bind(pins, p)
    verify(pins)
    result = dict(status='verified_original_domain_native_completion_api_terminal_unavailable',
        checked_utc=datetime.now(timezone.utc).isoformat(), **checked, unit=docs['launch']['unit'],
        invocation_id=e['invocation_id'], wrapper=e['wrapper'], child=e['child'],
        validation_sha256=sha(NAMES['producer']), original_wall_seconds=e['wall_seconds'],
        producer_declared_source_files=len(v['source_hashes']), producer_source_archive=str(archive),
        producer_source_archive_sha256=sha(archive), all_bulk_sources_freshly_rehashed_here=False,
        original_producer_rehashed_sources_before_terminal=True, original_journal=str(journal),
        source_hashes=pins, independent_atom_readback_complete=False, scientific_eligibility=False,
        gpu=False, new_predictions=0,
        scope='Exact original native wrapper/child/configuration, complete original start/end payloads '
              'and invocation-linked manager records prove native zero completion. Original API64749 '
              'is unavailable and its terminal remains null, never imputed zero. Original producer '
              'bulk-source hashes are retained as expected digests; a separate full reader must '
              'freshly replay shards/CIF atoms. No extraction restart or biological acceptance.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'wrapper', 'child']}, indent=2))


if __name__ == '__main__':
    main()
