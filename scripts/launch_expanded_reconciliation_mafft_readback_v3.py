#!/usr/bin/env python3
"""Launch a bounded independent MAFFT-reconciliation readback after exact systemd evidence validation.

This recovery is for user-systemd services that record invocation-linked start and
terminal resource messages but do not expose a _PID/_CMDLINE record for the
service's executable.  It never treats terminal unit fields alone as proof.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil

from audit_busco_gene_copies import ROOT, sha


PRODUCER_LAUNCH = Path('metadata/expanded_reconciliation_mafft_recovery_20261006_v2_launch.json')
PRODUCER_RECEIPT = Path('results/orthology/expanded-reconciliation-mafft-recovery-20261006-v2/receipt.json')
READER_PLAN = Path('metadata/expanded_reconciliation_mafft_recovery_20261006_v2_readback_plan.json')
FAILED_WAITER = Path('metadata/expanded_reconciliation_mafft_recovery_20261006_v2_readback_wait_plan.json')
TRANSPORT = Path('metadata/expanded_reconciliation_mafft_recovery_20261006_v3_readback_dependency_transport.json')
LAUNCH = Path('metadata/expanded_reconciliation_mafft_recovery_20261006_v3_readback_launch.json')
UNIT = 'fungal-expanded-reconciliation-mafft-recovery-readback-20261006-v3.service'
OUTPUT = Path('results/orthology/expanded-reconciliation-mafft-recovery-readback-20261006-v2')


def write_new(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def unit_state(unit):
    fields = ('ActiveState', 'Result', 'ExecMainStatus')
    result = {}
    for field in fields:
        result[field] = subprocess.check_output(
            ['systemctl', '--user', 'show', unit, '-p', field, '--value'], text=True).strip()
    return result


def journal_rows(unit):
    raw = subprocess.check_output(['journalctl', '--user', '-u', unit, '-o', 'json', '--no-pager'], text=True)
    return raw, [json.loads(line) for line in raw.splitlines() if line]


def verify_producer():
    launch = json.loads(PRODUCER_LAUNCH.read_text())
    expected = 'Started ' + launch['unit'] + ' - ' + ' '.join(launch['cmdline']) + '.'
    state = unit_state(launch['unit'])
    if state != {'ActiveState': 'inactive', 'Result': 'success', 'ExecMainStatus': '0'}:
        raise RuntimeError('Producer terminal state is not successful: ' + repr(state))
    raw, rows = journal_rows(launch['unit'])
    own = [row for row in rows if row.get('USER_INVOCATION_ID') == launch['invocation_id']]
    starts = [row for row in own if row.get('MESSAGE') == expected]
    resources = [row for row in own if row.get('CPU_USAGE_NSEC') and 'Consumed ' in (row.get('MESSAGE') or '')]
    failures = [row for row in own if 'Main process exited' in (row.get('MESSAGE') or '') or 'Failed with result' in (row.get('MESSAGE') or '')]
    if len(starts) != 1 or len(resources) != 1 or failures:
        raise RuntimeError('Original invocation journal is incomplete or indicates failure')
    receipt = json.loads(PRODUCER_RECEIPT.read_text())
    if receipt.get('status') != 'complete_native_mafft_reconciliation_pending_independent_output_readback':
        raise RuntimeError('Completed native producer receipt is required')
    if receipt.get('statistics') != {'species': 526, 'genes': 5815847} or not receipt.get('source_inputs_unchanged'):
        raise RuntimeError('Unexpected producer receipt scope')
    result = ROOT / receipt['result']
    for name, digest in receipt['mandatory_artifacts'].items():
        if sha(result / name) != digest:
            raise RuntimeError('Changed mandatory producer artifact: ' + name)
    return launch, state, raw, starts[0], resources[0], receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    launch, state, raw, start, resource, receipt = verify_producer()
    report = {
        'status': 'passed_systemd_invocation_linked_mafft_producer_completion_validation',
        'checked_utc': datetime.now(timezone.utc).isoformat(),
        'producer_launch': str(PRODUCER_LAUNCH),
        'producer_launch_sha256': sha(PRODUCER_LAUNCH),
        'producer_receipt': str(PRODUCER_RECEIPT),
        'producer_receipt_sha256': sha(PRODUCER_RECEIPT),
        'terminal_state': state,
        'invocation_id': launch['invocation_id'],
        'exact_start_message': start['MESSAGE'],
        'completion_resource_message': resource['MESSAGE'],
        'completion_cpu_usage_nsec': resource['CPU_USAGE_NSEC'],
        'native_journal_sha256': hashlib.sha256(raw.encode()).hexdigest(),
        'failed_v2_waiter_preserved': str(FAILED_WAITER),
        'mandatory_artifacts': receipt['mandatory_artifacts'],
        'scope': 'Exact producer invocation ID, start command message, terminal resource message, successful unit state, immutable native receipt and mandatory-artifact hashes checked. The earlier waiter failure is preserved; no producer output is changed and no biological reconciliation conclusion is accepted.'
    }
    if args.validate_only:
        print(json.dumps(report, indent=2))
        return
    for path in (TRANSPORT, LAUNCH):
        if path.exists():
            raise FileExistsError(path)
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))) * 1024
    if available < 96 * 2**30 or shutil.disk_usage(ROOT).free < 1024 * 2**30:
        raise RuntimeError('Readback resource preflight failed')
    report['source_hashes'] = {str(path): sha(path) for path in (PRODUCER_LAUNCH, PRODUCER_RECEIPT, READER_PLAN, FAILED_WAITER, Path(__file__))}
    write_new(TRANSPORT, report)
    command = [str(ROOT / '.cache/envs/orthofinder/bin/python'), 'scripts/readback_expanded_reconciliation_mafft_recovery_v1.py', '--plan', str(READER_PLAN)]
    systemd = ['systemd-run', '--user', '--collect', '--unit=' + UNIT, '--working-directory=' + str(ROOT),
               '--property=CPUQuota=400%', '--property=MemoryMax=64G', '--property=MemorySwapMax=0',
               '--setenv=PYTHONUNBUFFERED=1', '--setenv=OPENBLAS_NUM_THREADS=1', '--setenv=OMP_NUM_THREADS=1',
               '--setenv=MKL_NUM_THREADS=1', *command]
    subprocess.run(systemd, check=True)
    pid = 0
    for _ in range(50):
        pid = int(subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'MainPID', '--value'], text=True))
        if pid:
            break
        time.sleep(0.1)
    if pid <= 0:
        raise RuntimeError('Readback service did not expose a PID')
    process = psutil.Process(pid)
    if process.cmdline() != command:
        raise RuntimeError('Readback service command differs')
    cgroup = subprocess.check_output(['systemctl', '--user', 'show', UNIT, '-p', 'ControlGroup', '--value'], text=True).strip()
    limits = {key: (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / key).read_text().strip()
              for key in ('cpu.max', 'memory.max', 'memory.swap.max')}
    expected_limits = {'cpu.max': '400000 100000', 'memory.max': str(64 * 2**30), 'memory.swap.max': '0'}
    if limits != expected_limits:
        raise RuntimeError('Unexpected readback cgroup limits: ' + repr(limits))
    write_new(LAUNCH, {'status': 'launched_bounded_independent_mafft_reconciliation_readback_after_explicit_systemd_validation',
                       'unit': UNIT, 'pid': pid, 'created': process.create_time(), 'cmdline': command,
                       'transport': str(TRANSPORT), 'transport_sha256': sha(TRANSPORT), 'reader_plan': str(READER_PLAN),
                       'reader_plan_sha256': sha(READER_PLAN), 'actual_cgroup_limits': limits,
                       'systemd_launch_command': systemd, 'scientific_eligibility': False,
                       'scope': 'Recovery launch for one independent native-output readback. No producer/source mutation, model fitting, structure prediction or biological inference.'})
    print(json.dumps({'unit': UNIT, 'pid': pid, 'launch': str(LAUNCH), 'transport': str(TRANSPORT)}, indent=2))


if __name__ == '__main__':
    main()
