#!/usr/bin/env python3
"""Record complete hashed stages and exact original completion journals; tolerate null structured journal messages."""
import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import psutil
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    historical = plan.get('historical_source_hashes', {})
    historical_seen = {}

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)

    verify()
    receipts = {}
    for label, spec in plan['evidence'].items():
        path = Path(spec['path'])
        if str(path) not in bindings:
            raise ValueError('Unpinned evidence: ' + str(path))
        record = json.loads(path.read_text())
        for field, value in spec['expected'].items():
            if record[field] != value:
                raise ValueError('Incomplete evidence: ' + label + ':' + field)
        for name, digest in record.get('artifacts', {}).items():
            artifact = str(path.parent / name)
            if sha(artifact) != digest:
                raise ValueError('Changed artifact: ' + artifact)
            bindings[artifact] = digest
        for field in ['source_hashes', 'source_pins', 'pins']:
            for source, digest in record.get(field, {}).items():
                if source in historical:
                    if historical[source] != digest:
                        raise ValueError('Historical source hash differs: ' + source)
                    historical_seen[source] = digest
                    continue
                if source in bindings and bindings[source] != digest:
                    raise ValueError('Inconsistent binding: ' + source)
                bindings[source] = digest
        receipts[label] = record
    for rule in plan['links']:
        if receipts[rule['from']][rule['field']] != bindings[rule['to']]:
            raise ValueError('Receipt lineage differs: ' + str(rule))

    services = []
    for path in plan['launches']:
        launch = json.loads(Path(path).read_text())
        unit = launch.get('unit', launch.get('service'))
        try:
            proc = psutil.Process(launch['pid'])
            if proc.create_time() == launch['created'] and proc.status() != psutil.STATUS_ZOMBIE:
                raise RuntimeError('Recorded producer remains live: ' + unit)
        except psutil.NoSuchProcess:
            pass
        command = ['systemctl', '--user', 'show', unit, '-p', 'LoadState',
                   '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus']
        state = dict(line.split('=', 1) for line in subprocess.check_output(command, text=True).splitlines())
        if {k: state[k] for k in ['ActiveState', 'Result', 'ExecMainStatus']} != {
                'ActiveState': 'inactive', 'Result': 'success', 'ExecMainStatus': '0'}:
            raise RuntimeError('Unsuccessful terminal observation: ' + unit)
        journal = subprocess.check_output(['journalctl', '--user', '-u', unit,
                                           '-o', 'json', '--no-pager'], text=True)
        entries = [json.loads(line) for line in journal.splitlines()]
        exact = [r for r in entries if r.get('_PID') == str(launch['pid'])
                 and r.get('_CMDLINE') == ' '.join(launch['cmdline'])]
        if not exact:
            raise ValueError('No journal evidence for captured process: ' + unit)
        invocations = {r['_SYSTEMD_INVOCATION_ID'] for r in exact}
        resources = [r for r in entries if r.get('USER_INVOCATION_ID') in invocations
                     and r.get('CPU_USAGE_NSEC')]
        if not resources or any('Main process exited' in (r.get('MESSAGE') or '')
                                or 'Failed with result' in (r.get('MESSAGE') or '') for r in entries):
            raise ValueError('Missing completion or failed invocation: ' + unit)
        services.append(dict(launch=path, observed_state=state,
                             journal_sha256=__import__('hashlib').sha256(journal.encode()).hexdigest(),
                             captured_process_messages=len(exact),
                             completion_resource_records=[{
                                 k: r.get(k) for k in ['__REALTIME_TIMESTAMP', 'CPU_USAGE_NSEC',
                                                     'MEMORY_PEAK', 'MEMORY_SWAP_PEAK']
                             } for r in resources],
                             terminal_evidence='Captured process is absent; exact PID/command journal completion and complete hashed artifacts/proofs checked. '
                             'LoadState=not-found means a collected transient unit: its default success/0 fields alone do not prove exit status.'))
    verify()
    output = Path(plan['output'])
    if set(historical_seen) != set(historical):
        raise ValueError('Historical source hash was not encountered')
    result = dict(status=plan['completed_status'], checked_utc=datetime.now(timezone.utc).isoformat(),
                  summary=plan['summary'], scientific_eligibility=False,
                  services=services, source_hashes=bindings,
                  historical_source_hashes=historical_seen, scope=plan['scope'])
    with output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], output=str(output),
                          bound_sources=len(bindings), services_checked=len(services))))


if __name__ == '__main__':
    main()
