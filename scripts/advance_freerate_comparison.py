#!/usr/bin/env python3
"""Wait for full FreeRate optimization, audit it, then compare selected estimates."""
import argparse
import fcntl
import json
import subprocess
import sys
import time
from pathlib import Path
from audit_busco_gene_copies import ROOT, sha
from advance_paired_path_uncertainty import process_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    config_hash = sha(args.config)
    output = ROOT / config['controller_output']
    output.mkdir(parents=True, exist_ok=True)
    lock = (output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (output / 'receipt.json').exists():
        raise FileExistsError('Controller already completed')

    def verify():
        if sha(args.config) != config_hash:
            raise ValueError('Controller configuration changed')
        for path, digest in config['pinned_files'].items():
            if sha(ROOT / path) != digest:
                raise ValueError('Pinned file changed: ' + path)

    verify()
    while True:
        current = process_identity(config['producer_pid'])
        if current is None or current['start_ticks'] != config['producer_start_ticks'] or current['state'] == 'Z':
            break
        print('waiting_for_full_optimization', config['producer_pid'], flush=True)
        time.sleep(30)
    verify()
    producer_path = ROOT / config['optimization'] / 'receipt.json'
    producer = json.loads(producer_path.read_text())
    if (producer['status'] != 'complete_flagged_freerate_diagnostic_execution'
            or producer['source_fits'] != config['expected_source_fits']
            or producer['diagnostic_fits'] != 4 * config['expected_source_fits']):
        raise ValueError('Full optimization completion not established')
    audit_command = [sys.executable, 'scripts/audit_freerate_optimization.py',
                     '--diagnostics', config['optimization'], '--free', config['free'],
                     '--output', config['optimization_audit']]
    with (output / 'audit.log').open('w') as log:
        subprocess.run(audit_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    audit_path = ROOT / config['optimization_audit'] / 'receipt.json'
    audit = json.loads(audit_path.read_text())
    if (audit['status'] != 'passed_full_freerate_optimization_diagnostics'
            or audit['source_fits'] != config['expected_source_fits']
            or audit['diagnostic_receipt_sha256'] != sha(producer_path)):
        raise ValueError('Full optimization audit not established')
    command = [sys.executable, 'scripts/compare_site_rate_heterogeneity.py']
    for name in ('gamma', 'gamma_audit', 'free', 'free_audit', 'inputs', 'optimization', 'optimization_audit'):
        command.extend(['--' + name.replace('_', '-'), config[name]])
    command.extend(['--output', config['comparison_output']])
    with (output / 'comparison.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    result_path = ROOT / config['comparison_output'] / 'receipt.json'
    result = json.loads(result_path.read_text())
    if result['status'] != 'complete_matched_rate_heterogeneity_comparison' or result['fits'] != config['expected_source_fits']:
        raise ValueError('Full selected-fit comparison not established')
    receipt = {'status': 'complete_freerate_comparison_handoff', 'config_sha256': config_hash,
               'script_sha256': sha(Path(__file__)), 'producer_receipt_sha256': sha(producer_path),
               'audit_receipt_sha256': sha(audit_path), 'comparison_receipt_sha256': sha(result_path),
               'commands': [audit_command, command],
               'interpretation': 'Execution and audit completed; independent selected-fit comparison readback remains a separate step.'}
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
