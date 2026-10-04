#!/usr/bin/env python3
"""Qualify actual GDB native-write capture and clean normal inferior exit."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha, write_json
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists(); a.output.mkdir(exist_ok=False)
    pins = {}; modules = project_sources(pins, [Path(__file__),
        Path('scripts/run_retained_factor_watchpoint_debugger_v2.py'),
        Path('scripts/retained_factor_watchpoint_control_v2.py')])
    for name in [sys.executable, '/usr/bin/gdb', '/usr/bin/prlimit']: bind(pins, name)
    outcomes = {}
    for case in ['native_write', 'normal_exit']:
        folder = a.output/case; folder.mkdir()
        target = folder/'target'; debugger = folder/'debugger'; metadata = target/'watchpoints.json'
        target_receipt = target/'receipt.json'; plan_path = folder/'plan.json'
        plan = dict(output=str(target), debugger_output=str(debugger),
            watchpoint_metadata=str(metadata), target_receipt=str(target_receipt),
            target_expected_status='completed_artificial_watchpoint_normal_exit_target',
            target_command=[sys.executable, 'scripts/retained_factor_watchpoint_control_v2.py',
                '--metadata', str(metadata), '--receipt', str(target_receipt), '--case', case],
            artificial_control=True, pins=dict(pins),
            scope='Actual native GDB control on artificial words only; no original factors or biological inference.')
        write_json(plan_path, plan); receipt_path = folder/'driver_receipt.json'
        with (folder/'stdout.log').open('x') as out, (folder/'stderr.log').open('x') as err:
            proc = subprocess.run([sys.executable, 'scripts/run_retained_factor_watchpoint_debugger_v2.py',
                '--plan', str(plan_path), '--receipt', str(receipt_path)], stdout=out, stderr=err, timeout=60)
        assert proc.returncode == 0
        receipt = json.loads(receipt_path.read_text()); verify(receipt['source_hashes'])
        assert receipt['debugger_exit_code'] == 0 and len(receipt['installed_watchpoints']) == 2
        if case == 'native_write':
            assert receipt['status'] == 'captured_native_watchpoint_positive_control'
            assert receipt['native_write_hit'] and receipt['stopped_program_counters']
            assert not receipt['target_receipt_present'] and not receipt['normal_exit_inspection_skipped']
            text = (debugger/'gdb_stdout.log').read_text()
            assert '__memset' in text and 'Old value = ' in text and 'New value = 0' in text
        else:
            assert receipt['status'] == 'completed_artificial_watchpoint_normal_exit_control'
            assert not receipt['native_write_hit'] and not receipt['stopped_program_counters']
            assert receipt['target_receipt_present'] and receipt['normal_exit_inspection_skipped']
            assert 'No registers.' not in (debugger/'gdb_stderr.log').read_text()
            assert json.loads(target_receipt.read_text())['words_preserved']
        outcomes[case] = dict(status=receipt['status'], debugger_exit_code=0,
            driver_exit_code=0, installed_watchpoints=2,
            native_write_hit=receipt['native_write_hit'],
            normal_exit_inspection_skipped=receipt['normal_exit_inspection_skipped'],
            receipt=str(receipt_path), receipt_sha256=sha(receipt_path))
    for path in a.output.rglob('*'):
        if path.is_file(): bind(pins, path)
    verify(pins)
    result = dict(status='passed_actual_native_watchpoint_debugger_v2_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), outcomes=outcomes,
        actual_native_debugger_controls=2, artificial_control=True,
        transitive_project_source_modules=len(modules), source_hashes=pins,
        original_factor_diagnostic_runs=0, original_jobs_restarted=False,
        biological_fits=0, scientific_eligibility=False, gpu=False,
        scope='Two actual GDB controls: capture deliberately written artificial word and skip live-register '
              'inspection after normal inferior exit. No original factors, math/source repair, original '
              'native retry, unprotected corruption attribution or biological acceptance.')
    write_json(a.receipt, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
