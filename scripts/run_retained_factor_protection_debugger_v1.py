#!/usr/bin/env python3
"""Retain a fresh protected-factor diagnostic and native debugger evidence."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha, write_json
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); verify(plan['pins'])
    pins = dict(plan['pins']); bind(pins, args.plan)
    root = Path(plan['debugger_output']); root.mkdir(exist_ok=False)
    target_receipt = root / 'target_receipt.json'
    expressions = ['set pagination off', 'set confirm off', 'set debuginfod enabled off',
        'set disable-randomization off', 'set auto-load python-scripts off',
        'handle SIGSEGV stop print nopass', 'run',
        'printf "\\nPROJECT_FAULT_ADDRESS=%p\\n", $_siginfo._sifields._sigfault.si_addr',
        'thread apply all bt 30', 'info registers', 'x/8i $pc', 'info sharedlibrary']
    command = ['/usr/bin/gdb', '-nx', '--batch']
    for expression in expressions: command += ['-ex', expression]
    command += ['--args', sys.executable, 'scripts/diagnose_retained_timing_protected_factors_v7.py',
                '--plan', str(args.plan), '--receipt', str(target_receipt)]
    write_json(root/'debugger_command.json', command)
    with (root/'gdb_stdout.log').open('x') as out, (root/'gdb_stderr.log').open('x') as err:
        child = subprocess.Popen(command, stdout=out, stderr=err, env=dict(os.environ, DEBUGINFOD_URLS=''))
        process = psutil.Process(child.pid)
        identity = dict(pid=process.pid, created=process.create_time(), cmdline=process.cmdline())
        write_json(root/'debugger_process.json', identity)
        code = child.wait()
    target_root = Path(plan['output'])
    text = (root/'gdb_stdout.log').read_text()
    fault = 'received signal SIGSEGV' in text
    addresses = re.findall(r'PROJECT_FAULT_ADDRESS=(0x[0-9a-fA-F]+)', text)
    ranges_path = target_root/'factor_protection_ranges.json'
    ranges = json.loads(ranges_path.read_text()) if ranges_path.exists() else {}
    protected_faults = [tree for tree, interval in ranges.items()
                       if any(interval['start'] <= int(address, 16) < interval['end_exclusive']
                              for address in addresses)]
    if target_receipt.exists():
        assert not fault
        target = json.loads(target_receipt.read_text())
        assert target['status'] == 'completed_protected_factor_timing_probe_diagnostic_v7'
        assert target['protected_factor_interior_pages']
        status = ('completed_protected_factor_probe_without_native_fault' if code == 0
                  else 'completed_target_probe_debugger_post_exit_commands_require_review')
    elif fault and protected_faults:
        status = 'captured_native_fault_in_protected_source_factor'
    else:
        status = 'retained_debugger_or_inferior_failure_requires_review'
    for folder in [root, target_root]:
        for path in folder.rglob('*'):
            if path.is_file(): bind(pins, path)
    verify(pins)
    result = dict(status=status, checked_utc=datetime.now(timezone.utc).isoformat(),
        plan=str(args.plan), plan_sha256=sha(args.plan), debugger_exit_code=code,
        debugger_identity=identity, target_receipt_present=target_receipt.exists(),
        native_sigsegv_observed=fault, fault_addresses=addresses,
        protected_factor_faults=protected_faults, source_hashes=pins,
        original_jobs_restarted=False, production_tolerance_changed=False,
        scientific_eligibility=False, biological_fits=0, gpu=False,
        scope='Separate process-local protected-factor probe under GDB. Raw native backtrace, '
              'instruction, fault address, protection ranges and exact debugger/inferior identities '
              'retained. Protected fault is diagnostic evidence, not a repair or attribution of '
              'the original unprotected timing failure. Partial edge pages remain writable. '
              'No original restart, tolerance relaxation, optimizer, MCMC or biological acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ['source_hashes', 'debugger_identity']}, indent=2))


if __name__ == '__main__':
    main()
