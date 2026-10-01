#!/usr/bin/env python3
"""Revalidate original completion identities without treating JSON key order as content.

Checks the closed archive's actual hash and original PID/create/command journal
messages and completion resource records. Full scientific inputs are not
rehashed here; their prior complete closure remains a separate requirement.
Raw journal serialization differences are diagnostics, never silently rewritten.
"""
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def audit(completion, journal_output):
    completion = Path(completion); closed = json.loads(completion.read_text())
    archive_path = Path(closed['full_hash_archive'])
    assert sha(archive_path) == closed['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    assert len(archive['source_hashes']) == closed['bound_source_hashes']
    assert len(archive['services']) == closed['exact_process_journals_checked']
    assert all(closed.get(k) == v for k, v in archive['summary'].items())
    bindings = {str(completion): sha(completion), str(archive_path): sha(archive_path)}
    for key in ['producer_receipt', 'independent_readback', 'source_plan']:
        path = closed[key]; assert sha(path) == closed[key + '_sha256']; bindings[path] = sha(path)
    output = Path(journal_output); output.mkdir(exist_ok=False, parents=True)
    services = []
    for old in archive['services']:
        lp = Path(old['launch']); launch = json.loads(lp.read_text())
        assert sha(lp) == archive['source_hashes'][str(lp)]; bindings[str(lp)] = sha(lp)
        unit = launch.get('unit', launch.get('service'))
        try:
            process = psutil.Process(launch['pid'])
            assert process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE, unit
        except psutil.NoSuchProcess:
            pass
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), unit
        raw = subprocess.check_output(['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True)
        entries = [json.loads(line) for line in raw.splitlines()]
        exact = [r for r in entries if r.get('_PID') == str(launch['pid']) and r.get('_CMDLINE') == ' '.join(launch['cmdline'])]
        assert exact and len(exact) == old['captured_process_messages'], unit
        invocations = {r['_SYSTEMD_INVOCATION_ID'] for r in exact}
        resources = [r for r in entries if r.get('USER_INVOCATION_ID') in invocations and r.get('CPU_USAGE_NSEC')]
        fields = ['__REALTIME_TIMESTAMP', 'CPU_USAGE_NSEC', 'MEMORY_PEAK', 'MEMORY_SWAP_PEAK']
        assert resources and [{k: r.get(k) for k in fields} for r in resources] == old['completion_resource_records'], unit
        original = [r for r in entries if r.get('_SYSTEMD_INVOCATION_ID') in invocations or r.get('USER_INVOCATION_ID') in invocations]
        assert not any('Main process exited' in r.get('MESSAGE', '') or 'Failed with result' in r.get('MESSAGE', '') for r in original), unit
        path = output / (lp.stem + '.jsonl'); path.write_text(raw); bindings[str(path)] = sha(path)
        canonical = '\n'.join(json.dumps(r, sort_keys=True, separators=(',', ':')) for r in entries)
        services.append(dict(launch=str(lp), unit=unit, pid=launch['pid'], created=launch['created'], cmdline=launch['cmdline'],
            original_invocations=sorted(invocations), observed_state=state, captured_process_messages=len(exact),
            completion_resource_records=old['completion_resource_records'], original_completion_fields_equal=True,
            current_journal=str(path), current_journal_sha256=sha(path),
            current_canonical_journal_sha256=hashlib.sha256(canonical.encode()).hexdigest(),
            closure_raw_journal_sha256=old['journal_sha256'], raw_serialization_equal=sha(path) == old['journal_sha256']))
    for path, digest in bindings.items(): assert sha(path) == digest, path
    return dict(status='passed_closed_original_journal_semantic_revalidation', checked_utc=datetime.now(timezone.utc).isoformat(),
        completion=str(completion), completion_sha256=sha(completion), closed_archive=str(archive_path), closed_archive_sha256=sha(archive_path),
        exact_process_journals_checked=len(services), raw_serialization_differences=sum(not r['raw_serialization_equal'] for r in services),
        services=services, source_hashes=bindings, scientific_eligibility=False, scope=__doc__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completion', type=Path, required=True); parser.add_argument('--journal-output', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    result = audit(args.completion, args.journal_output)
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'exact_process_journals_checked', 'raw_serialization_differences']}))
